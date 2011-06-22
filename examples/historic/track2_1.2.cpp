#include "lp_lib.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "ANN.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include <sys/stat.h>

using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::ostringstream;
using boost::shared_array;
using boost::shared_ptr;

const int VERSION = 1; // Version of program
const int SUBVERSION = 2; // Subversion of program
const string DATE_STRING = "2010/06/15"; // Date of current version
const string AUTHOR_STRING = "Frederik Kaster"; // Author(s) of current version

#define isnan( x ) ( (x) != (x) )

// C/C++ implementation of the module track2.m (see comments in this file)
//
// Track cells in a greedy fashion (time slice by time slice) by 
// minimizing a simple energy function (mainly based on distances) via
// Integer Linear Programming.
//
// Implementation is geared towards trivial (data-parallel) parallelization.
//
//  Requires the following 3rd party libraries:
// - LP_SOLVE for mixed-integer programming 
//   http://lpsolve.sourceforge.net/5.5/
// - HDF5 for input and output data storage
//   http://www.hdfgroup.org/HDF5/
// - ANN for efficient nearest-neighbor search
//   http://www.cs.umd.edu/~mount/ANN/
//
// Frederik Orlando Kaster, 2010 <frederik.kaster@iwr.uni-heidelberg.de>
//

// Parameters for the matching process 
class MatchingPars {
public:
  // number of nearest neighbors to be considered as children
  int knn; 
  // maximum distance between a cell and its potential children
  double maxDist; 
  // cost for a split
  double cSplit; 
  
  double cDeath; // cost for a cell death
  
  double cAppear; // cost for an appearance

  double xyzRatio; // Ratio of lateral resolution to axial resolution (typically > 1)
  
  unsigned long flagFeatures; // flag determining which features shall be used
  
  MatchingPars(); // Constructor with default arguments
  
  int readOptions(int nOpts, char* opts[]); // read parameters from string array
};

MatchingPars::MatchingPars() : knn(6),
                               maxDist(25),
                               cSplit(100),
                               cDeath(1000),
                               cAppear(10000),
			       xyzRatio(1.),
                               flagFeatures(0x200) {
}

// Read parameters from the input options: see function printHelp() for
// explanations
int MatchingPars::readOptions(int nOpts, char* opts[]){
  if( nOpts % 2 ){
    cerr << "Options must be passed as key-value pairs" << endl;
    return 1;
  }
  nOpts/=2;
  const string allowedKeys("kmsdafz");
  for( int iOpt=0; iOpt<nOpts; ++iOpt ){
    string key(opts[2*iOpt]);
    string value(opts[2*iOpt+1]);
    size_t optNr = 0;
    if( key.size()!=2 || key[0]!='-' || 
        (optNr=allowedKeys.find_first_of(key[1]))==string::npos ){
      cerr << "Invalid key for option no. " << iOpt << endl;
      return 1;
    } else {
      std::istringstream valueReader(value);
      switch(optNr){
      case 0 : valueReader >> knn; break;
      case 1 : valueReader >> maxDist; break;
      case 2 : valueReader >> cSplit; break;
      case 3 : valueReader >> cDeath; break;
      case 4 : valueReader >> cAppear; break;
      case 5 : valueReader >> std::hex >> flagFeatures >> std::dec; break;
      case 6 : valueReader >> xyzRatio; break;
      }
      if( valueReader.fail() ){
        cerr << "Invalid value for option no. " << iOpt << endl;
        return 1;
      }
    }
  }
  return 0;
} 

// Check whether the file exists
int fileExists( const std::string& fileName ) {
  int output=0;
  if( !fileName.empty() ){
    struct stat buf;
    output = (stat( fileName.c_str(), &buf ) == 0);
  }
  return output;
}

// Print the usage instructions for this program 
void printHelp() {
  cout << "TRACK2 : Simple tracking of presegmented cells\n"
       << "Version " << VERSION << "." << SUBVERSION << ", " << DATE_STRING
       << ", " << AUTHOR_STRING << "\n\n"
       << "Usage: track2 FILESLICE1 FILESLICE2 FILEOUT [OPTIONS]\n"
       << "FILESLICE1/2 : files with cell coordinates for the two time "
          "slices (HDF5)\n"
       << "FILEOUT : file where the matching info is written (HDF5)\n"
       << "Note : if FILEOUT exists, the program does nothing\n"
       << "(existing output files are not overwritten)\n"
       << "Available options (examples):\n"
       << "-k 3 : set the number of nearest neighbors to 3 (default: 6)\n"
       << "-m 15 : set the maximum parent/child distance to 15 (default: "
          "25)\n"
       << "-s 30 : set the cost for a split to 30 (default: 100)\n"
       << "-d 200 : set the cost for a cell death to 200 (default: 1000)\n"
       << "-a 1000 : set the cost of an appearance to 1000 (default: "
          "10000)\n"
       << "-f 10300 : set the flags of which features and use for tracking to "
          "0x10300 (default: 0x00200)\n\n"
       << "Supported flags:\n" 
       << "00100 : center of binary mask\n"
       << "00200 : center of mass (i.e. weighted mask)\n"
       << "10000 : total intensity of cell\n"
       << "-z 10 : set the ratio of axial (z) to lateral (xy) resolution to "
          "10 (default: 1)\n" 
       << endl; 
}

// General interface class for reading a data file with the features 
// required for the tracking
//
class FileReader {
public:

// Construct from file name and Z / XY resolution ratio
  FileReader( const std::string& fileName, double xyzRatio=1. );

// Read the input features from the encapsulated data file. 
// Upon successful completion, nCells contains the number of cells,
// nFeats contains the number of features and feats contains the nFeats*nCells
// feature data (first all features for cell no. 1 stacked behind each other,
// then all features for cell no. 2, ...).
  int readInputs(int& nCells, int& nFeats, 
		 boost::shared_array<double>& feats);
  
// Read-only access to file name
  std::string getFileName() const;

// Get a vector of all valid cell labels for this file
  void getValidLabels(std::vector<int>& validLabels) const;

// Get the Z / XY resolution ratio
  double getXYZRatio() const;
protected:
  virtual int readInputsSpecific(int& nCells, int& nFeats,
                                 boost::shared_array<double>& feats) = 0;
// Default implementation: returns all numbers between 1 and nCells
  virtual void getValidLabelsSpecific(std::vector<int>& validLabels) const;
private:
  std::string fileName_;
  int nCells_;
// Z to XY resolution ratio (should be >= 1.)
  double xyzRatio_; 
};

FileReader::FileReader( const string& fileName, double xyzRatio ) : fileName_(fileName),
  nCells_(0), xyzRatio_(xyzRatio) {
}

std::string FileReader::getFileName() const {
  return fileName_;
}

double FileReader::getXYZRatio( ) const {
  return xyzRatio_;
}

void FileReader::getValidLabels( std::vector<int>& validLabels ) const {
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


// Read the files from the EMBL Philipp Keller datafiles: see
// http://www.embl.de/digitalembryo/

class FileReaderEmbl : public FileReader {
public:
// Constructor
  FileReaderEmbl( const std::string& fileName, double xyzRatio=1. );

protected:   
// Read the input coordinates from a specified HDF5 file in the 
// EMBL format (see http://www.embl.de/digitalembryo/)
   virtual int readInputsSpecific(int& nCells, int& nFeats,
                                 boost::shared_array<double>& feats);
 
};

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
  if( dataDims[0]!=nFeats ){
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

// Reader for the input coordinates from the HCI format for digital embryo
// files (see 
// https://gorgonzola.iwr.uni-heidelberg.de/intern/wiki/index.php/Interface )
class FileReaderHci : public FileReader {
public:
// Constructor with initialization of admissible named features
  FileReaderHci( const std::string& fileName, double xyzRatio=1. );
// Determine which features shall be used for the tracking
  int setFlag( unsigned long flag );
protected:
  virtual int readInputsSpecific( int& nCells, int& nFeats, 
                                  shared_array<double>& feats );
  virtual void getValidLabelsSpecific( std::vector<int>& validLabels ) const;
private:
  // Named Feature structure that contains all necessary information for
  // retrieving a feature and the corresponding information from an HDF5 file
  struct NamedFeature {
    std::string name_; // e.g. "position"
    int nSubFeats_;    // how many subfeatures this feature comprises, e.g.
                       // 3 for a position
    shared_array<double> weights_; // weights to be used for this feature
  };

  // Named flags that can be used to decide which data to extract from
  // the file
  // key (long) -- Flag used in the "featureconfig" data element and showing
  //               which features have been computed
  // value (NamedFeature) -- Corresponding feature to be stored             
  std::map<unsigned long, NamedFeature> namedFlags_;
    
  unsigned long flag_; // flag deciding which features shall be extracted

  // Vector containing 1's for valid cell labels and 0's for invalid labels
  std::vector<char> labelValidity_;
};

FileReaderHci::FileReaderHci( const string& fileName, double xyzRatio ) : 
  FileReader(fileName, xyzRatio), namedFlags_(), flag_(0), labelValidity_() {
  NamedFeature position = {"position", 3, shared_array<double>(new double[3]) };
  position.weights_[0] = 1.;
  position.weights_[1] = 1.;
  position.weights_[2] = xyzRatio;
  namedFlags_[0x100] = position;
  NamedFeature com = {"com", 3, shared_array<double>(new double[3]) };
  com.weights_[0] = 1.;
  com.weights_[1] = 1.;
  com.weights_[2] = xyzRatio;
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
  const string nLabelsName = "supervoxels";
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
  if( (err<0) || (dataType!=H5T_INTEGER) || (dims1d!=nLabels) ){
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
            (dims1d < currFeat.nSubFeats_) ){
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

// Compute the optimum matching by optimizing a linear integer program.
//
// The inputs are usually acquired by the readInputs() function and the
// MatchingPars::readOptions() function.
// Upon successful completion, nMoves / nSplits / nDeaths / nAppears
// denote the numbers of motion / mitosis / death or disappearance /
// appearance events.
//
// moves : nMoves * 2 matrix with cell indices in the first and second
//         slice
// splits : nSplits * 3 matrix, with the cell indices of the mother cells in
//          the first row, and the indices of the daughter cells in the 
//          second and third row
// deaths : nDeaths * 1 matrix, with the indices of the disappearing cells
//          from the first time slice
// appears : nAppears * 1 matrix, with the indices of the appearing cells 
//           from the second time slice
//
// All matrices are vectorized in row-major order (as usual in C).
// solveFlag shows the output flag of the LP_SOLVE operation, see
//
// http://lpsolve.sourceforge.net/5.5/solve.htm
//
// for the interpretation of the various integral return values.
//
// totCost holds the total costs of the optimum matching.


int computeMatch( int nCells1, shared_array<double> feats1, 
                  const vector<int>& validLabels1,
                  int nCells2, shared_array<double> feats2,
                  const vector<int>& validLabels2,
                  int nFeats,
                  vector<int>& moves, vector<int>& splits,
                  vector<int>& deaths, vector<int>& appears,
                  double& totCost, int& solveFlag,
                  const MatchingPars& pars ){
  // to make sure, clear the inputs
  moves.clear();
  splits.clear();
  deaths.clear();
  appears.clear();
  lprec* lp = NULL;
  ANNpointArray pts2 = NULL;
  ANNpoint currQ = NULL;
  ANNidxArray currNn = NULL;
  ANNdistArray currDists = NULL;
  int output = 0;
  try{
    int output = 1;
    // initialize the linear program -> it will have nCells1+nCells2 linear
    // constraint and an (unknown) number of variables 
    // we start with an empty program and add the columns of the constraint
    // matrix successively
    int nConstr = nCells1 + nCells2;
    int nVars = 0; // will be incremented : number of variables
    lp = make_lp(nConstr, nVars);
    if( lp==NULL ){
      throw "Could not initialize linear program";
    }
    // construct a KD tree for finding the nearest neighbors among the
    // potential daughters and decreasing the number of variables in the
    // linear program
    const double eps = 0; // allowed error for kNN searches
    pts2 = annAllocPts( nCells2, nFeats );
    if( pts2==NULL ){
      throw "Allocation of points for kd-tree failed";
    }
    for( int i=0; i<nCells2; ++i ){
      ANNpoint cp = pts2[i];
      for( int j=0; j<nFeats; ++j ){
        cp[j] = feats2[i*nFeats+j];
      }
    }
    currDists = new ANNdist[pars.knn];
    shared_ptr<ANNkd_tree> kdTree( new ANNkd_tree( pts2, nCells2, nFeats ) );
    currQ = annAllocPt(nFeats);
    currNn = new ANNidx[pars.knn];
    if( currDists==NULL || kdTree==NULL || currQ==NULL || currNn==NULL ){
      throw "Allocation failure";
    }
    // Storage for building the sparse representation of the constraint
    // matrix
    REAL moveColData[3] = {0.,0.,0.};
    REAL splitColData[4] = {0.,0.,0.,0.};
    int moveRowNrs[3] = {0,0,0};
    int splitRowNrs[4] = {0,0,0,0};
    // store the representation of move and split moves as 3 subsequent
    // integers : 
    // MOTHER_INDEX DAUGHTER_INDEX -1 (for a move)
    // MOTHER_INDEX DAUGHTER1_INDEX DAUGHTER2_INDEX (for a split)
    vector<int> variables;
    for( int i=0; i<nCells1; ++i ){
      for( int j=0; j<nFeats; ++j ){
        currQ[j]=feats1[i*nFeats+j];
      } 
      int totPts = kdTree->annkFRSearch( currQ, pars.maxDist*pars.maxDist,
	 pars.knn, currNn, currDists, eps );
      totPts = (totPts < pars.knn) ? totPts : pars.knn;
      // add a column for every potential move
      for( int k=0; k<totPts; ++k ){
        if( isnan( currDists[k] ) ){
          continue;
        }
        moveColData[0] = currDists[k]-pars.cDeath-pars.cAppear;
        moveColData[1] = 1.;
        moveColData[2] = 1.;
        moveRowNrs[0] = 0;
        moveRowNrs[1] = i+1;
        moveRowNrs[2] = nCells1+currNn[k]+1;
        unsigned char out = add_columnex( lp, 3, moveColData, moveRowNrs );
        if( !out ){
          throw "Problem while adding move column to LP";
        }
        nVars++;
        out = set_binary( lp, nVars, TRUE );
        variables.push_back( i );
        variables.push_back( currNn[k] );
        variables.push_back( -1 );
        if( !out ){
          throw "Problem while constraining variable to be binary";
        }
      }
      // add a constraint column for every potential split
      for( int k=0; k<totPts; ++k ){
        if( isnan(currDists[k]) ){
          continue;
        }
        for( int l=k+1; l<totPts; ++l ){
          if( isnan(currDists[l]) ){
            continue;
          }
          splitColData[0] = currDists[k]+currDists[l]+pars.cSplit
				-pars.cDeath-2*pars.cAppear;
          splitColData[1] = 1.;
          splitColData[2] = 1.;
          splitColData[3] = 1.;
          splitRowNrs[0] = 0;
          splitRowNrs[1] = i+1;
          splitRowNrs[2] = nCells1+currNn[k]+1;
          splitRowNrs[3] = nCells1+currNn[l]+1;
          unsigned char out = add_columnex( lp, 4, splitColData, splitRowNrs );
          if( !out ){
            throw "Problem while adding split column to LP";
          }
          nVars++;
          out = set_binary( lp, nVars, TRUE );
          variables.push_back( i );
          variables.push_back( currNn[k] );
          variables.push_back( currNn[l] );
          if( !out ){
            throw "Problem while constraining variable to be binary";
          }
        }
      }
    }
    // set the right-hand side vector of the constraints
    shared_array<REAL> rhs( new REAL[nConstr+1] );
    rhs[0] = 0; // will be ignored
    unsigned char out = 1;
    for( int i=1; i<=nCells1; ++i ){
      rhs[i] = 1;
      out = set_constr_type( lp, i, LE );
      if( !out ){
        break;
      }
    }
    if( !out ){
      throw "Error while setting constraint type";
    }
    for( int i=nCells1+1; i<=nConstr; ++i ){
      rhs[i]=1;
      out = set_constr_type( lp, i, LE );
      if(!out){
        break;
      }
    }
    if( !out ){
      throw "Error while setting constraint type";
    }
    set_rh_vec( lp, rhs.get() );
    set_minim( lp );
    solveFlag = solve( lp );
    switch(solveFlag){
    case NOMEMORY: throw "Out of memory while solving LP program"; break;
    case OPTIMAL: cout << "Proper convergence" << endl; break;
    case SUBOPTIMAL: cout << "Suboptimal solution found" << endl; break;
    case INFEASIBLE: throw "Model is infeasible"; break;
    case UNBOUNDED: throw "Model is unbounded"; break;
    case DEGENERATE: cout << "Model is degenerate" << endl; break;
    case NUMFAILURE: throw "Numerical failure during LP"; break;
    case USERABORT: throw "LP aborted by user"; break;
    case TIMEOUT: throw "Aborted because of time-out"; break;
    case PRESOLVED: cout << "Model can be presolved perfectly" << endl; break;
    case PROCFAIL: throw "Branch and bound failed"; break;
    case PROCBREAK: throw "Break at first / break at value"; break;
    case FEASFOUND: cout << "Feasible branch and bound solution found" << endl;
                    break;
    case NOFEASFOUND: throw "No feasible branch and bound solution found"; 
                      break;
    }
    totCost = get_working_objective( lp ) + pars.cDeath * nCells1
              + pars.cAppear * nCells2;
    shared_array<REAL> finalVars( new REAL[nVars] );
    out = get_variables( lp, finalVars.get() );
    if( !out ){
      throw "Error while extracting the final variables";
    }
    shared_array<bool> accountCells1(new bool[nCells1]);
    for( int i=0; i<nCells1; ++i ){
      accountCells1[i]=false;
    }
    shared_array<bool> accountCells2(new bool[nCells2]);
    for( int i=0; i<nCells2; ++i ){
      accountCells2[i]=false;
    }
    for( int i=0; i<nVars; ++i ){
      if( finalVars[i] == 0 ){
        continue;
      }
      if( variables[3*i+2] == -1 ){ // move
        moves.push_back(validLabels1[variables[3*i]]);
        moves.push_back(validLabels2[variables[3*i+1]]);
        accountCells1[variables[3*i]]=true;
        accountCells2[variables[3*i+1]]=true;
      } else { // split
        splits.push_back(validLabels1[variables[3*i]]);
        splits.push_back(validLabels2[variables[3*i+1]]);
        splits.push_back(validLabels2[variables[3*i+2]]);
        accountCells1[variables[3*i]]=true;
        accountCells2[variables[3*i+1]]=true;
        accountCells2[variables[3*i+2]]=true;
      }
    }
    // the cells unaccounted for must be deaths or appearances
    for( int i=0; i<nCells1; ++i ){
      if( !accountCells1[i] ){
        deaths.push_back( validLabels1[i] );
      }
    }
    for( int i=0; i<nCells2; ++i ){
      if( !accountCells2[i] ){
        appears.push_back( validLabels2[i] );
      }
    }
  } catch ( char* str ){
    cerr << str << endl;
    output=1; 
  }
  if( currQ ){
    annDeallocPt(currQ);
  }
  if( pts2 ){
    annDeallocPts(pts2);
  }
  if( lp ){
    delete_lp(lp);
  }
  if( currNn ){
    delete[] currNn;
  }
  if( currDists ){
    delete[] currDists;
  }
  annClose();
  return output;
}

/* Write the outputs of the linear program (moves and splits, cell deaths,
   cell appearances, the total cost and the return flag of the linear 
   program) to an HDF5 file 
*/

int writeOutputs( const char* fileName, const MatchingPars& pars,
                  const vector<int>& moves, const vector<int>& splits, 
                  const vector<int>& deaths, const vector<int>& appears,
                  double totCost, int solveFlag ){
  hid_t fid = H5Fcreate( fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
  if( fid<0 ){
    cerr << "Creation of output file " << fileName << " failed" << endl;
    return 1;
  }
  if( !moves.empty() ){
    const string movesName = "/Moves";
    hsize_t dimsMoves[2];
    dimsMoves[0] = moves.size()/2;
    dimsMoves[1] = 2;
    shared_array<int> movesBuffer( new int[moves.size()] );
    for( unsigned i=0; i<moves.size(); ++i ){
      movesBuffer[i] = moves[i];
    }
    herr_t err = H5LTmake_dataset_int( fid, movesName.c_str(), 2, dimsMoves,
                                       movesBuffer.get() );
    if( err<0 ){
      cerr << "Dataset " << movesName << " could not be written" << endl;
      return 1;
    }
  }
  if( !splits.empty() ){
    const string splitsName = "/Splits";
    hsize_t dimsSplits[2];
    dimsSplits[0] = splits.size()/3;
    dimsSplits[1] = 3;
    shared_array<int> splitsBuffer( new int[splits.size()] );
    for( unsigned i=0; i<splits.size(); ++i ){
      splitsBuffer[i] = splits[i];
    }
    herr_t err = H5LTmake_dataset_int( fid, splitsName.c_str(), 2, dimsSplits, 
      splitsBuffer.get() );
    if( err<0 ){
      cerr << "Dataset " << splitsName << " could not be written" << endl;
      return 1;
    }
  }
  if( !deaths.empty() ){
    const string deathsName = "/Disappearances";
    hsize_t dimDeaths = deaths.size();
    shared_array<int> deathsBuffer( new int[deaths.size()] );
    for( unsigned i=0; i<deaths.size(); ++i ){
      deathsBuffer[i] = deaths[i];
    }
    herr_t err = H5LTmake_dataset_int( fid, deathsName.c_str(), 1, &dimDeaths, 
      deathsBuffer.get() );
    if( err<0 ){
      cerr << "Dataset " << deathsName << " could not be written" << endl;
      return 1;
    }
  }
  if( !appears.empty() ){
    const string appearsName = "/Appearances";
    hsize_t dimAppears = appears.size();
    shared_array<int> appearsBuffer( new int[appears.size()] );
    for( unsigned i = 0; i<appears.size(); ++i ){
      appearsBuffer[i] = appears[i];
    }
   herr_t err = H5LTmake_dataset_int(fid, appearsName.c_str(), 1, &dimAppears, 
      appearsBuffer.get() );
    if ( err<0 ){
      cerr << "Dataset " << appearsName << " could not be written" << endl;
      return 1;
    }
  }
  const string versionName = "/Version";
  hsize_t dimScalar = 1;
  herr_t err = H5LTmake_dataset_int(fid, versionName.c_str(), 1, &dimScalar,
     &VERSION);
  if( err < 0 ){
    cerr << "Version tag could not be written" << endl;
    return 1;
  }
  const string subversionName = "/Subversion";
  err = H5LTmake_dataset_int(fid, subversionName.c_str(), 1, &dimScalar,
     &SUBVERSION);
  if( err < 0 ){
    cerr << "Subversion tag could not be written" << endl;
    return 1;
  } 
  const string solveFlagName = "/Output flag of LP_SOLVE";
  err = H5LTmake_dataset_int(fid, solveFlagName.c_str(), 1, &dimScalar, 
    &solveFlag );
  if( err<0 ){
    cerr << "Dataset " << solveFlagName << " could not be written" << endl;
    return 1;
  }
  const string totCostName = "/Total costs";
  err = H5LTmake_dataset_double(fid, totCostName.c_str(), 1, &dimScalar, 
    &totCost );
  if( err<0 ){
    cerr << "Dataset " << totCostName << " could not be written" << endl;
    return 1;
  }
  // Write the matching parameters
  const string knnName = "/Number of nearest neighbors as possible children";
  herr_t err1 = H5LTmake_dataset_int(fid, knnName.c_str(), 1, &dimScalar, 
                                     &pars.knn);
  const string maxDistName = "/Maximum parent-child distance";
  herr_t err2 = H5LTmake_dataset_double(fid, maxDistName.c_str(), 1, &dimScalar,
                                        &pars.maxDist);
  const string cSplitName = "/Cost of cell split";
  herr_t err3 = H5LTmake_dataset_double(fid, cSplitName.c_str(), 1, &dimScalar,
                                        &pars.cSplit);
  const string cDeathName = "/Cost of cell death";
  herr_t err4 = H5LTmake_dataset_double(fid, cDeathName.c_str(), 1, &dimScalar,
                                        &pars.cDeath);
  const string cAppearName = "/Cost of cell appearance";
  herr_t err5 = H5LTmake_dataset_double(fid, cAppearName.c_str(), 1, &dimScalar,
                                        &pars.cAppear);
  const string cFlagsName = "/Feature extraction flags";
  herr_t err6 = H5LTmake_dataset(fid, cFlagsName.c_str(), 1, &dimScalar,
                                 H5T_NATIVE_ULONG, &pars.flagFeatures);
  const string xyzRatioName = "/Ratio between lateral (XY) and axial (Z) resolution";
  herr_t err7 = H5LTmake_dataset_double(fid, xyzRatioName.c_str(), 1, &dimScalar,
                                        &pars.xyzRatio);
  bool errorOutput= (err1<0) || (err2<0) || (err3<0) || (err4<0) || (err5<0) ||
                    (err6<0) || (err7<0);
  if( errorOutput ){
    cerr << "Error while writing matching parameters" << endl;
  }
  err = H5Fclose( fid );
  if( err<0 ){
    cerr << "Error while closing file" << endl;
  }
  return 0;
}

// Run the program
int main(int argc, char *argv[]){
  if( argc<4 ){
    printHelp();
    return 1;
  }
  if( fileExists( argv[3] ) ){
    cout << "Output file " << argv[3] << " already exists: abort" << endl;
    return 0;
  }
  MatchingPars pars;
  int output = pars.readOptions( argc - 4, argv + 4 );
  if( output ){
    cerr << "Options wrongly specified\n";
    printHelp();
    return 2;
  }
  int nCells1 = 0, nCells2 = 0;
  int nFeats1 = 0, nFeats2 = 0;
  shared_array<double> feats1, feats2;
  FileReaderHci reader1( argv[1] );
  reader1.setFlag( pars.flagFeatures );
  output = reader1.readInputs( nCells1, nFeats1, feats1 );
  if( output ){
    cerr << "Error reading coordinates from " << argv[1] << endl;
    return 5;
  }
  vector<int> validLabels1;
  reader1.getValidLabels( validLabels1 );
  FileReaderHci reader2( argv[2] );
  reader2.setFlag( pars.flagFeatures );
  output = reader2.readInputs( nCells2, nFeats2, feats2 );
  if( output ){
    cerr << "Error reading coordinates from " << argv[2] << endl;
    return 6;
  }
  vector<int> validLabels2;
  reader2.getValidLabels( validLabels2 );
  vector<int> moves, deaths, splits, appears;
  double totCost = 0;
  int solveFlag = 0;
  cout << "Computation of matches between " << argv[1] << " and " 
       << argv[2] << " started" << endl;
  output = computeMatch( nCells1, feats1, validLabels1, nCells2, feats2, 
                         validLabels2, nFeats1, moves, splits,
                         deaths, appears, totCost, solveFlag, pars );
  if( output ){
    cerr << "Error during main computation" << endl;
    return 7;
  }
  output = writeOutputs( argv[3], pars, moves, splits, deaths, appears,
                         totCost, solveFlag );
  if( output ){
    cerr << "Error while writing results to " << argv[3] << endl;
    return 8;
  }
  cout << "Computation ended successfully" << endl;
  return 0;
}


