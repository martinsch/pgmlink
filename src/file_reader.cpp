#include "file_reader.h"
#include "common.h"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "hdf5.h"
#include "hdf5_hl.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::ostringstream;
using boost::shared_array;

namespace Tracking {

FileReader::FileReader( const string& fileName, double xyzRatio ) : fileName_(fileName),
  nCells_(0), xyzRatio_(xyzRatio) {
}

string FileReader::getFileName() const {
  return fileName_;
}

double FileReader::getXYZRatio( ) const {
  return xyzRatio_;
}

void FileReader::getValidLabels( vector<int>& validLabels ) const {
  validLabels.clear();
  getValidLabelsSpecific( validLabels );
}

void FileReader::getValidLabelsSpecific(vector<int>& validLabels) const {
  validLabels.resize( nCells_ );
  for( int i=0; i<nCells_; ++i ){
    validLabels[i] = i + 1;
  }
}

int FileReader::readInputs(int& nCells, int& nFeats,
                           shared_array<double>& feats){
  int output=1;
  if( fileExists( fileName_ ) ){
    output=readInputsSpecific( nCells, nFeats, feats);
  }
  if( output ){
    nCells=nFeats=0;
    feats.reset();
  }
  nCells_ = nCells;
  return output;
}

FileReaderEmbl::FileReaderEmbl( const string& fileName, double xyzRatio ) : FileReader(fileName, xyzRatio){
}

int FileReaderEmbl::readInputsSpecific( int& nCells, int& nFeats, 
                                        shared_array<double>& feats ){
  nFeats = 3;
  const string datasetName = "/COM Coordinates";
  string fileName = getFileName();
  hid_t fid = H5Fopen( fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
  if( fid<0 ){
    cerr << "File " << fileName << " cannot be opened" << endl;
    return 1;
  }
  int dataRank = 0;
  herr_t err = H5LTget_dataset_ndims( fid, datasetName.c_str(), &dataRank );
  if( err<0 ){
    cerr << "Dataset '" + datasetName + "' does not exist in " << fileName << endl;
    return 1;
  }
  if( dataRank!=2 ){
    cerr << "Dataset '" + datasetName + "' in " << fileName << " has wrong rank"
         << endl;
    return 1;
  }
  hsize_t dataDims[2] = { 0, 0 };
  H5T_class_t dataType = H5T_FLOAT;
  size_t dataTypeSize = 0;
  err = H5LTget_dataset_info( fid, datasetName.c_str(), dataDims, &dataType,
                              &dataTypeSize );
  if( err<0 ){
    cerr << "Error while getting properties of '" + datasetName + "' in file "
         << fileName << endl;
    return 1;
  }
  if( dataDims[0]!=static_cast<unsigned>(nFeats) ){
    cerr << "Data set " + datasetName + " in file " + fileName + " of wrong "
      "size" << endl;
    return 1;
  }
  if( dataType!=H5T_FLOAT ){
    cerr << "Data set " + datasetName + " in file " + fileName + " of wrong type"
   << endl;
    return 1;
  }
  nCells = dataDims[1];
  long nData=nFeats*nCells;
  feats.reset( new double[nData] );
  shared_array<double> featsBuf( new double[nData] ); // wrong order, rows 
                                                      // instead of cols
  err = H5LTread_dataset_double( fid, datasetName.c_str(), featsBuf.get() );
  shared_array<double> featsRatio( new double[nFeats] );
  featsRatio[nFeats-1] = getXYZRatio();
  for( int i = nFeats-2; i>=0; --i ){
    featsRatio[i] = 1.;
  }
  long idx=0;
  for( int iC=0; iC<nCells; iC++ ){
    for( int iF=0; iF<nFeats; iF++ ){
      feats[idx++] = featsBuf[iF*nCells+iC] * featsRatio[iF];
    }
  } 
  featsBuf.reset();
  if( err<0 ){
    cerr << "Error reading data from "+datasetName+" in file "+fileName << endl;
    return 1;
  }
  err = H5Fclose( fid );
  if( err < 0 ){
    cerr << "Error while closing file" << endl;
    return 1;
  }
  return 0;
}

FileReaderHci::FileReaderHci( const string& fileName, double xyzRatio ) : 
  FileReader(fileName, xyzRatio), namedFlags_(), flag_(0), labelValidity_() {
  NamedFeature position = {"position", 3, shared_array<double>(new double[3]) };
  position.weights_[0] = xyzRatio;
  position.weights_[1] = 1.;
  position.weights_[2] = 1.;
  namedFlags_[0x100] = position;
  NamedFeature com = {"com", 3, shared_array<double>(new double[3]) };
  com.weights_[0] = xyzRatio;
  com.weights_[1] = 1.;
  com.weights_[2] = 1.;
  namedFlags_[0x200] = com;
  NamedFeature intensity = {"intensity", 1, shared_array<double>(new double[1])};
  intensity.weights_[0] = 1.;
  namedFlags_[0x10000] = intensity;
}

int FileReaderHci::setFlag( unsigned long flag ) {
  int output = 0;
  unsigned long currFlag=1; // used for extracting bits of flag from the right
  // check whether all bits in flag correspond to an entry of namedFlags_
  unsigned long copyFlag = flag;
  while( copyFlag ){
    if( (currFlag & copyFlag)!=0 ){
      if( namedFlags_.find(currFlag)==namedFlags_.end() ){ // not found
        output = 1;
        break;
      }
      copyFlag -= currFlag;
    }
    currFlag <<= 1;
  }
  if( output==0 ){
    flag_ = flag;
  } else {
    flag_ = 0;
  }
  return output;
}

int FileReaderHci::readInputsSpecific( int& nCells, int& nFeats, 
                                       shared_array<double>& feats ){
  // find out which cells are valid
  labelValidity_.clear();
  string fileName = getFileName();
  hid_t fid = H5Fopen( fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
  if( fid < 0 ){
    cerr << "File " << fileName << " cannot be opened" << endl;
    return 1;
  }
  int dataRank = 0;
  const string groupName = "/features/";
  // find out how many labels there are in the file  
  const string nLabelsName = "labelcount";
  string datasetName = groupName + nLabelsName;
  herr_t err = H5LTget_dataset_ndims( fid, datasetName.c_str(),
                                      &dataRank );
  if( (err<0) || (dataRank!=1) ){
    cerr << "Error while getting rank for " << datasetName << " from "
         << fileName << endl;
    return 1;
  }
  hsize_t dims1d = 0;
  H5T_class_t dataType = H5T_INTEGER;
  size_t dataTypeSize = 0;
  err = H5LTget_dataset_info( fid, datasetName.c_str(), &dims1d, &dataType,
                              &dataTypeSize );
  if( (err<0) || (dataType!=H5T_INTEGER) || (dims1d!=1) ){
    cerr << "Error while getting info for " << datasetName << " from "
         << fileName << endl;
    return 1;
  }
  int nLabels = 0;
  err = H5LTread_dataset_int( fid, datasetName.c_str(), &nLabels );
  if( err<0 ){
    cerr << "Error while reading " << datasetName << " from " << fileName 
         << endl;
    return 1;
  }
  // find which labels are valid
  vector<char> newLabelValidity( nLabels );
  const string validLabelsName = "labelcontent";
  datasetName = groupName + validLabelsName;
  err = H5LTget_dataset_ndims( fid, datasetName.c_str(), &dataRank );
  if( (err<0) || (dataRank!=1) ){
    cerr << "Error while getting rank for " << datasetName << " from "
         << fileName << endl;
    return 1;
  }
  err = H5LTget_dataset_info( fid, datasetName.c_str(), &dims1d, &dataType,
                              &dataTypeSize );
  if( (err<0) || (dataType!=H5T_INTEGER) || (dims1d!=static_cast<unsigned>(nLabels)) ){
    cerr << "Error while getting info for " << datasetName << " from " 
         << fileName << endl;
    return 1;
  }
  shared_array<unsigned short> shortBuffer(new unsigned short[nLabels]);
  err = H5LTread_dataset( fid, datasetName.c_str(), H5T_NATIVE_USHORT,
                          shortBuffer.get() );
  if( err<0 ){
    cerr << "Error while reading " << datasetName << " from " << fileName 
         << endl;
    return 1;
  }
  nCells = 0;
  for( int i=0; i<nLabels; ++i ){
    newLabelValidity[i] = (char) shortBuffer[i];
    if( newLabelValidity[i] ){
      nCells++;
    }
  }
  // now that the number of cells is known, find the number of features based
  // on the flags
  unsigned long currFlag = 1;
  unsigned long copyFlag = flag_;
  vector<unsigned long> allFlags;
  while( copyFlag ){
    if( (currFlag & copyFlag)!=0 ){
      allFlags.push_back( currFlag ); 
      copyFlag -= currFlag;
    }
    currFlag <<= 1; 
  }
  nFeats = 0;
  int nFlags = allFlags.size();
  for( int i = 0; i < nFlags; ++i ){
    nFeats += namedFlags_[ allFlags[i] ].nSubFeats_;
  }
  feats.reset( new double[nFeats*nCells] );
  const string featuresString = "/features/";
  ostringstream oss;
  long arrayIdx = 0; // current index of the feature array 
  for( int i = 0; i < nLabels; ++i ){
    if( newLabelValidity[i] ){
      oss.str("");
      oss << i + 1 << "/"; 
      for( vector<unsigned long>::const_iterator it=allFlags.begin(); 
           it!=allFlags.end(); ++it ){
        const NamedFeature& currFeat = namedFlags_[ *it ];
        datasetName = featuresString + oss.str() + currFeat.name_;
        err = H5LTget_dataset_ndims( fid, datasetName.c_str(), &dataRank );
        if( (err < 0) || (dataRank != 1) ){
          cerr << "Error while getting rank for " << datasetName << " from "
               << fileName << endl;
          return 1;
        }  
        err = H5LTget_dataset_info( fid, datasetName.c_str(), &dims1d,
                                    &dataType, &dataTypeSize );
        if( (err < 0) || (dataType != H5T_FLOAT) || 
            (dims1d < static_cast<unsigned>(currFeat.nSubFeats_)) ){
          cerr << "Error while getting info for " << datasetName << " from "
               << fileName << endl;
          return 1;
        }
        shared_array<double> doubleBuffer( new double[dims1d] );
        err = H5LTread_dataset( fid, datasetName.c_str(), H5T_NATIVE_DOUBLE,
                                doubleBuffer.get() );
        if( err < 0 ){
          cerr << "Error while reading " << datasetName << " from " << fileName
               << endl;
          return 1;
        }
        for( int j = 0; j < currFeat.nSubFeats_; ++j ){
          feats[arrayIdx++] = doubleBuffer[j] * currFeat.weights_[j];
        }
      }
    }
  }
  err = H5Fclose( fid );
  if( err < 0 ){
    cerr << "Error while closing file " + fileName << endl;
    return 1;
  }
  labelValidity_.swap( newLabelValidity );
  return 0;
}

void FileReaderHci::getValidLabelsSpecific( vector<int>& validLabels ) const {
  int potLabelSize = labelValidity_.size();
  validLabels.reserve( potLabelSize );
  for( int i=0; i<potLabelSize; ++i ){
    if( labelValidity_[i] ){
      validLabels.push_back(i + 1);
    }
  }
}



}
