#ifdef USE_CPLEX
#include <ilcplex/cplex.h>
#else
#include "lp_lib.h"
#endif /* USE_CPLEX */

#include "hdf5.h"
#include "hdf5_hl.h"
#include "ANN.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>
#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include <sys/stat.h>
#include <ctime>

using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::ostringstream;
using boost::shared_array;
using boost::shared_ptr;

const int VERSION = 1; // Version of program
const int SUBVERSION = 3; // Subversion of program
const string DATE_STRING = "2010/06/28"; // Date of current version
const string AUTHOR_STRING = "Frederik Kaster"; // Author(s) of current version

enum ReturnCode { SUCCESSFUL, READING_ERROR, COMPUTATION_ERROR, WRITING_ERROR };

const size_t chunkSize = 50; // For creating a chunked data set

#define isnan( x ) ( (x) != (x) )

#ifdef USE_CPLEX
const string solverName = "CPLEX 12.1";
#else
const string solverName = "LP_SOLVE 5.5";
#endif /* USE_CPLEX */
// the following size constant must be large enough to hold all strings we may want to write to an HDF5 file
const size_t stringBufferSize = 80;

// C/C++ implementation of the module track2.m (see comments in this file)
//
// Track cells in a greedy fashion (time slice by time slice) by 
// minimizing a simple energy function (mainly based on distances) via
// Integer Linear Programming.
//
// Implementation is geared towards trivial (data-parallel) parallelization.
//
//  Requires the following 3rd party libraries:
// - LP_SOLVE or CPLEX for mixed-integer programming 
//   http://lpsolve.sourceforge.net/5.5/
//   http://www.ilog.com/products/cplex/
// - HDF5 for input and output data storage
//   http://www.hdfgroup.org/HDF5/
// - ANN for efficient nearest-neighbor search
//   http://www.cs.umd.edu/~mount/ANN/
//
// Frederik Orlando Kaster, 2010 <frederik.kaster@iwr.uni-heidelberg.de>
//
// If CPLEX shall be used as optimization engine, this program has to be built
// with the USE_CPLEX flag.
//
// Changelog:
// 1.3 -- optional support for using CPLEX as a MIP engine
//        OpenMPI support

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
       << "Usage: track2 TASKFILE [OPTIONS]\n"
       << "TASKFILE : text file with the names of all files on which TRACK2 shall"
          " operate\n"
       << "These files are assumed to be in subsequent order\n"
       << "Upon successful completion, the tracking information describing the\n"
       << "assocation between each pair of subsequent files will be written into the\n"
       << "second file of this pair, if necessary overwriting any previous tracking\n"
       << "information. Hence read-write access to the files listed in TASKFILE\n"
       << "must be possible\n"
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

class IlpSolver {

public:
  // Common interface of the ILP solvers:
  //
  // We assume that all variables must be binary, and that all
  // inequalities are of the <= form.
  // nVars - number of variables (columns of constraint matrix)
  // nConstr - number of constraints (rows of constraint matrix)
  // nNonZero - number of non-zero constraint matrix entries
  // costs - costs of setting a variable to 1 (length nVars)
  // rhs - right-hand sides of inequalities (length nConstr)
  // matbeg - index of the beginning of every column in the 
  //          coefficient array matval (ascending, length nVars)
  // matcnt - number of nonzero elements in each column
  // matind - row numbers of the coefficients in matval (length nNonZero)
  // matval - array of nonzero coefficients (length nNonZero)
  // finalCost - costs of the final IP solution
  // finalVars - final values of the variables at the optimum
  // solverOutput - solver-specific output flag
  // outputString - string explaining this output flag
  // 
  // The output determines whether the operation was successful (0), or
  // not (1)
  virtual int solve(  int nVars, int nConstr, int nNonZero,
                      const vector<double>& costs,
                      const vector<double>& rhs,
                      const vector<int>& matbeg,
                      const vector<int>& matcnt,
                      const vector<int>& matind,
                      const vector<double>& matval,
                      const string& problemString,
                      double& finalCost,
                      shared_array<double>& finalVars,
                      int& solverOutput,
                      string& outputString ) = 0;
};

#ifdef USE_CPLEX
class CplexSolver : public IlpSolver {
public:
   CplexSolver();
   ~CplexSolver();
   virtual int solve( int nVars, int nConstr, int nNonZero,
                      const vector<double>& costs,
                      const vector<double>& rhs,
                      const vector<int>& matbeg,
                      const vector<int>& matcnt,
                      const vector<int>& matind,
                      const vector<double>& matval,
                      const string& problemString,
                      double& finalCost,
                      shared_array<double>& finalVars,
                      int& solverOutput,
                      string& outputString );
private:
  CPXENVptr env_;
  // report a CPLEX error message to stderr
  void reportCplexError( int errcode );
  // deactivate copying
  CplexSolver(const CplexSolver&);
  CplexSolver& operator=(const CplexSolver&);
};

CplexSolver::CplexSolver() : env_(NULL) {
  int status;
  env_ = CPXopenCPLEX(&status);
  if(env_==NULL){
    reportCplexError( status );
    throw "Error while initializing CPLEX environment";
  }
}

CplexSolver::~CplexSolver(){
  int status = CPXcloseCPLEX(&env_);
  if( status ){
    fprintf(stderr, "Error while releasing CPLEX environment\n");
  }
}

void CplexSolver::reportCplexError(int errcode){
  char errbuffer[4096];
  CPXCCHARptr outptr = CPXgeterrorstring(env_, errcode, errbuffer);
  if( outptr!=NULL ){
    fprintf(stderr, "%s\n", errbuffer);
  } else {
    fprintf(stderr, "Unknown CPLEX error flag: %d\n", errcode); 
  }
}

// Solve an integer linear program using the CPLEX solver
int CplexSolver::solve( int nVars, int nConstr, int nNonZero,
                        const vector<double>& costs,
                        const vector<double>& rhs,
                        const vector<int>& matbeg,
                        const vector<int>& matcnt,
                        const vector<int>& matind,
                        const vector<double>& matval,
                        const string& problemString,
                        double& finalCost,
                        shared_array<double>& finalVars,
                        int& solverOutput,
                        string& outputString){
  int output = 1;
  CPXLPptr lp = NULL;
  try{
    int status = 0;
    lp = CPXcreateprob(env_, &status, problemString.c_str() );
    if( lp==NULL ){
      reportCplexError(status);
      throw "Error while creating LP problem " + problemString; 
    }
    shared_array<char> sense(new char[nConstr]);
    for( int c=0; c<nConstr; ++c ){
      sense[c] = 'L'; // <= constraint
    }
    shared_array<double> lb(new double[nVars]);
    shared_array<double> ub(new double[nVars]);
    for( int v=0; v<nVars; ++v ){
      lb[v] = 0.;
      ub[v] = 1.;
    }
    status = CPXcopylp(env_, lp, nVars, nConstr, CPX_MIN, &(costs[0]), &(rhs[0]),
                       sense.get(), &(matbeg[0]), &(matcnt[0]), &(matind[0]), 
                       &(matval[0]), lb.get(), ub.get(), NULL);
    if( status ){
      reportCplexError(status);
      throw "Error while filling ILP problem " + problemString;
    }
    shared_array<char> ctypes(new char[nVars]);
    for( int v=0; v<nVars; ++v ){
      ctypes[v] = CPX_BINARY;
    }
    status = CPXcopyctype(env_, lp, ctypes.get());
    if( status ){
      reportCplexError(status);
      throw "Error while setting variable types for problem " + problemString;
    }
    status = CPXmipopt(env_, lp);
    if( status ){
      reportCplexError(status);
      throw "Error while optimizing ILP problem " + problemString;
    }
    finalVars.reset( new double[nVars] );
    status = CPXsolution(env_, lp, &solverOutput, &finalCost, finalVars.get(), NULL,
      NULL, NULL);
    if( status ){
      reportCplexError(status);
      throw "Error while extracting the solution from the ILP problem " + problemString;
    }
    switch( solverOutput ){
    case CPXMIP_OPTIMAL: outputString = "Optimal solution found"; output = 0; break;
    case CPXMIP_OPTIMAL_INFEAS: outputString = "Optimal with unscaled infeasibilities";
      output = 0; break;
    case CPXMIP_OPTIMAL_TOL: outputString = "Optimal solution within tolerance"; 
      output = 0; break;
    case CPXMIP_TIME_LIM_FEAS: outputString = "Time limit exceeded, problem feasible";
      break;
    case CPXMIP_TIME_LIM_INFEAS: 
      outputString = "Time limit exceeded, no feasible solution"; break;
    case CPXMIP_NODE_LIM_FEAS: outputString = "Node limit exceeded, problem feasible";
      break;
    case CPXMIP_NODE_LIM_INFEAS: 
      outputString = "Node limit exceeded, no feasible solution"; break;
    case CPXMIP_UNBOUNDED: outputString = "Problem is unbounded"; break;
    case CPXMIP_SOL_LIM: outputString = "Limit on MIP solutions reached"; break;
    case CPXMIP_INFEASIBLE: outputString = "Problem is integer infeasible"; break;
    case CPXMIP_INForUNBD: outputString = "Problem infeasible or unbounded"; break;
    case CPXMIP_MEM_LIM_FEAS: outputString = "Memory limit exceeded, problem feasible";
      break;
    case CPXMIP_MEM_LIM_INFEAS: outputString = "Memory limit exceeded, no feasible "
      "solution"; break;
    case CPXMIP_ABORT_FEAS: outputString = "Aborted, problem feasible"; break;
    case CPXMIP_ABORT_INFEAS: outputString = "Aborted, no feasible solution"; break;
    case CPXMIP_FAIL_FEAS: outputString = "Failure, problem feasible"; break;
    case CPXMIP_FAIL_FEAS_NO_TREE: outputString = "Out of memory, no tree available, "
      "problem feasible"; break;
    case CPXMIP_FAIL_INFEAS: outputString = "Failure, no feasible solution"; break;
    case CPXMIP_FAIL_INFEAS_NO_TREE: outputString = "Out of memory, no tree available,"
      " no feasible solution"; break;
    default: outputString = "Unexpected error flag, consult CPLEX manual"; break;
    }
  } catch( string& str ){
    cerr << str << endl;
    output = 1;
  } 
  if( lp ){
    int out = CPXfreeprob( env_, &lp );
    if( out ){
      fprintf(stderr, "Error while releasing LP problem");
    }
  }
  return output;
}

#else /* USE_CPLEX */

class LpSolveSolver : public IlpSolver {
public:  
  virtual int solve(  int nVars, int nConstr, int nNonZero,
                      const vector<double>& costs,
                      const vector<double>& rhs,
                      const vector<int>& matbeg,
                      const vector<int>& matcnt,
                      const vector<int>& matind,
                      const vector<double>& matval,
                      const string& problemString,
                      double& finalCost,
                      shared_array<double>& finalVars,
                      int& solverOutput,
                      string& outputString );
};

// Solve an integer linear program using the LP_SOLVE solver
int LpSolveSolver::solve( int nVars, int nConstr, int nNonZero,
                          const vector<double>& costs,
                          const vector<double>& rhs,
                          const vector<int>& matbeg,
                          const vector<int>& matcnt,
                          const vector<int>& matind,
                          const vector<double>& matval,
                          const string& problemString,
                          double& finalCost,
                          shared_array<double>& finalVars,
                          int& solverOutput,
                          string& outputString ){
  int output = 1;
  lprec* lp = NULL;
  try {
    lp = make_lp(nConstr, nVars);
    if( lp == NULL ){
      throw "LP problem " + problemString + " cannot be initialized"; 
    }
    // fill the constraint matrix by columns
    for( int c=0; c<nVars; ++c ){
      int nValues = matcnt[c]+1;
      REAL *colVals = new REAL[nValues];
      int *rowNrs = new int[nValues];
      colVals[0] = costs[c];
      rowNrs[0] = 0;
      for( int v=1; v<nValues; ++v ){
        colVals[v] = matval[matbeg[c]+v-1];
        rowNrs[v] = matind[matbeg[c]+v-1]+1;
      }
      unsigned char success = set_columnex(lp, c+1, nValues, colVals, rowNrs);
      delete[] colVals;
      delete[] rowNrs;
      if( !success ){
        ostringstream errMsg;
        errMsg << "Error while setting column no. " << c  << " in LP problem " 
                  + problemString;
        throw errMsg.str().c_str();
      }
      success = set_binary( lp, c+1, TRUE );
      if( !success ){
        ostringstream errMsg;
        errMsg << "Error while setting column no. " << c << " in " + problemString +
                  " to be binary";
        throw errMsg.str();
      }
    }
    // set the right-hand side vector of the constraints
    shared_array<REAL> extendedRhs( new REAL[nConstr+1] );
    extendedRhs[0] = 0.;
    for( int r=1; r<=nConstr; ++r ){
      extendedRhs[r] = rhs[r-1];
      unsigned char success = set_constr_type( lp, r, LE );
      if( !success ){
        ostringstream out;
        out << "Error while setting constraint type in " + problemString 
               + " for row no. " << r;
        throw out.str();
      }
    }
    set_rh_vec(lp, extendedRhs.get() );
    set_minim( lp );
    string solveString;
    solverOutput = ::solve( lp );
    finalCost = get_working_objective( lp );
    finalVars.reset( new double[nVars] );
    unsigned char success = get_variables( lp, finalVars.get() );
    if( !success ){
      throw "Error while extracting the final variables from " + problemString;
    }
    switch( solverOutput ){
    case NOMEMORY: outputString = "Out of memory while solving LP program"; break;
    case OPTIMAL: outputString = "Proper convergence"; output = 0; break;
    case SUBOPTIMAL: outputString = "Suboptimal solution found"; output = 0; break;
    case INFEASIBLE: outputString =  "Model is infeasible"; break;
    case UNBOUNDED: outputString =  "Model is unbounded"; break;
    case DEGENERATE: outputString =  "Model is degenerate"; break;
    case NUMFAILURE: outputString =  "Numerical failure during LP"; break;
    case USERABORT: outputString =  "LP aborted by user"; break;
    case TIMEOUT: outputString =  "Aborted because of time-out"; break;
    case PRESOLVED: outputString =  "Model can be presolved perfectly"; break;
    case PROCFAIL: outputString =  "Branch and bound failed"; break;
    case PROCBREAK: outputString =  "Break at first / break at value"; break;
    case FEASFOUND: outputString =  "Feasible branch and bound solution found";
                    break;
    case NOFEASFOUND: outputString =  "No feasible branch and bound solution found"; 
                      break;
    default: outputString = "Unknown output code - consult LP_SOLVE manual"; break;
    }
  } catch ( string& str ){
    cerr << str << endl;
    output = 1;
  }
  if( lp ){
    delete_lp( lp );
  }
  return output;
}
#endif /* USE_CPLEX */

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
// appears : nAppears * 1 matrix, with the indices of the appearing cells          from the second time slice
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
                  int nFeats, const string& problemString, IlpSolver& solver,
                  vector<int>& moves, vector<int>& splits,
                  vector<int>& deaths, vector<int>& appears,
                  double& totCost, int& solveFlag, string& solverString,
                  const MatchingPars& pars){
  // to make sure, clear the inputs
  moves.clear();
  splits.clear();
  deaths.clear();
  appears.clear();
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
    // for explanation of the following variables, see the interface of IlpSolver::solve
    vector<double> costs;
    vector<double> rhs(nConstr, 1.);
    vector<int> matbeg;
    vector<int> matcnt;
    vector<int> matind;
    int nNonZero = 0;
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
      // add a column of the constraint matrix for every potential move
      for( int k=0; k<totPts; ++k ){
        if( isnan( currDists[k] ) ){
          continue;
        }
        double costValue =  currDists[k]-pars.cDeath-pars.cAppear ;
        costs.push_back( costValue );
        // the vector matval is filled later, as it contains only 1's
        matbeg.push_back(nNonZero);
        matcnt.push_back(2);
        matind.push_back(i);
        matind.push_back(nCells1+currNn[k]);
        nNonZero += 2;
        variables.push_back( i );
        variables.push_back( currNn[k] );
        variables.push_back( -1 );
      }
      // add a column of the constraint matrix for every potential split
      for( int k=0; k<totPts; ++k ){
        if( isnan( currDists[k] ) ){
          continue;
        }
        for( int l=k+1; l<totPts; ++l ){
          if( isnan( currDists[l] ) ){
            continue;
          }
          double costValue = currDists[k]+currDists[l]+pars.cSplit-pars.cDeath-2*pars.cAppear;
          costs.push_back( costValue );
          matbeg.push_back(nNonZero);
          matcnt.push_back(3);
          matind.push_back(i);
          matind.push_back(nCells1+currNn[k]);
          matind.push_back(nCells1+currNn[l]);
          nNonZero += 3;
          variables.push_back( i );
          variables.push_back( currNn[k] );
          variables.push_back( currNn[l] );
        }
      }
    }
    // set the values of the constraint matrix
    int nVars = matcnt.size();
    vector<double> matval(nNonZero, 1.);
    shared_array<double> finalVars;
    int solveOut = solver.solve(nVars, nConstr, nNonZero, costs, rhs, matbeg, matcnt,
      matind, matval, problemString, totCost, finalVars, solveFlag, solverString);  
    totCost += pars.cDeath * nCells1 + pars.cAppear * nCells2;
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
  if( currNn ){
    delete[] currNn;
  }
  if( currDists ){
    delete[] currDists;
  }
  annClose();
  return output;
}

template<typename T>
hid_t correspondingType(){
  cerr << "Attempt to use a type for which the corresponding HDF5 identifier has not"
          " been specified" << endl;
  exit(EXIT_FAILURE); /* Programming error = brute-force abort */
  return H5T_NATIVE_INT; /* just to avoid compiler warnings */
}

/* Specialize for the types that shall actually be used */
template<> hid_t correspondingType<int>() {
  return H5Tcopy(H5T_NATIVE_INT);
}

template<> hid_t correspondingType<double>(){
  return H5Tcopy(H5T_NATIVE_DOUBLE);
}

template<> hid_t correspondingType<unsigned long>(){
  return H5Tcopy(H5T_NATIVE_ULONG);
}

// make sure that all strings to be written to a file do not exceed the maximum
// length of stringBufferSize

template<> hid_t correspondingType<char>(){
  hid_t strType = H5Tcopy(H5T_C_S1);
  H5Tset_size( strType, stringBufferSize );
  H5Tset_strpad( strType, H5T_STR_NULLTERM );
  return strType;
}

/* Convenience function to either create or overwrite a new dataset in a file or
   group, depending on whether it already exists or not (for array datasets).
   loc_id - File or group ID of the parent group
   dsetName - Name of the data set
   nDims - rank = number of dimensions
   dims - extent along the different dimensions
   data - Vector to the raw data (must lie contiguously in memory)
   T - datatype of the data to be written, must be int or double
*/

template<typename T>
herr_t createOrOverwriteDataset( hid_t loc_id, const char* dsetName, int nDims,
                                 hsize_t *dims, const T* data ){
  hid_t type = correspondingType<T>();
  if( H5LTfind_dataset(loc_id, dsetName) ){ // overwrite an existing dataset
    hid_t did = H5Dopen2(loc_id, dsetName, H5P_DEFAULT);
    herr_t succ = H5Dset_extent(did, dims);
    if( succ<0 ){
      cerr << "Changing the size of dataset " << dsetName << " failed" << endl;
      return succ;
    }
    herr_t out = H5Dwrite(did, type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          data);
    H5Dclose(did);
    return out;
  } else { // create a new dataset
    hid_t props = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t *chunkDims = new hsize_t[nDims];
    chunkDims[0] = chunkSize;
    for( int i=1; i<nDims; ++i ){
      chunkDims[i] = dims[i];
    }
    herr_t out = H5Pset_chunk(props, nDims, chunkDims);
    hsize_t *maxDims = new hsize_t[nDims];
    for( int i=0; i<nDims; ++i ){
      maxDims[i] = H5S_UNLIMITED;
    }
    hid_t dspace = H5Screate_simple(nDims, dims, maxDims);
    hid_t did = H5Dcreate2(loc_id, dsetName, type, dspace, H5P_DEFAULT, props, H5P_DEFAULT);
    hid_t fileSpace = H5Dget_space( did );
    out = H5Dwrite(did, type, dspace, fileSpace, H5P_DEFAULT, data);
    H5Sclose(dspace);
    H5Dclose(did);
    delete[] chunkDims;
    delete[] maxDims;
    return out;
  }
}

/* Create or overwrite a scalar dataset
   loc_id - ID of file or group where the dataset shall reside
   dsetName - name of this dataset
   datum - content to be written
*/

template<typename T>
herr_t createOrOverwriteScalar( hid_t loc_id, const char* dsetName, const T* datum){
  hid_t type = correspondingType<T>();
  herr_t out = 0;
  if( H5LTfind_dataset(loc_id, dsetName) ){ // overwrite an existing dataset
    hid_t did = H5Dopen2(loc_id, dsetName, H5P_DEFAULT);
    out = H5Dwrite(did, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, datum);
    H5Dclose(did);
  } else { // create a new dataset
    hid_t dspace = H5Screate( H5S_SCALAR );
    hid_t did = H5Dcreate2(loc_id, dsetName, type, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    out = H5Dwrite(did, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, datum);
    H5Sclose(dspace);
    H5Dclose(did);
  }
  H5Tclose(type);
  return out;
}

/* Write the outputs of the linear program (moves and splits, cell deaths,
   cell appearances, the total cost and the return flag of the linear 
   program) to an existing HDF5 file 
*/

int writeOutputs( const char* fileName, const MatchingPars& pars,
                  const vector<int>& moves, const vector<int>& splits, 
                  const vector<int>& deaths, const vector<int>& appears,
                  double totCost, int solveFlag, const string& solveString,
                  double matchingTime ){
  hid_t fid = H5Fopen( fileName, H5F_ACC_RDWR, H5P_DEFAULT );
  if( fid<0 ){
    cerr << "Opening the output file " << fileName << " failed" << endl;
    return 1;
  }
  const string trackingGroupName = "tracking";
  hid_t gid; // ID of the "Tracking" group
  if( H5Lexists(fid, trackingGroupName.c_str(), H5P_DEFAULT) ){
    gid = H5Gopen2(fid, trackingGroupName.c_str(), H5P_DEFAULT);
  } else {
    gid = H5Gcreate2(fid, trackingGroupName.c_str(), 0, H5P_DEFAULT, H5P_DEFAULT);
  }
  if( gid<0 ){
    cerr << "Opening or creating the group " << trackingGroupName << " failed" << endl;
    return 1;
  }
  const string movesName = "Moves";
  if( !moves.empty() ){
    hsize_t dimsMoves[2];
    dimsMoves[0] = moves.size()/2;
    dimsMoves[1] = 2;
    herr_t err = createOrOverwriteDataset<int>( gid, movesName.c_str(), 2, dimsMoves,
                                                &(moves[0]) );
    if( err<0 ){
      cerr << "Dataset " << movesName << " could not be written" << endl;
      return 1;
    }
  } else {
    if( H5LTfind_dataset(gid, movesName.c_str()) ){
      herr_t status = H5Ldelete(gid, movesName.c_str(), H5P_DEFAULT);
      if( status<0 ){
        cerr << "Dataset " << movesName << " could not be written" << endl;
        return 1;
      }
    }
  }
  const string splitsName = "Splits";
  if( !splits.empty() ){
    hsize_t dimsSplits[2];
    dimsSplits[0] = splits.size()/3;
    dimsSplits[1] = 3;
    herr_t err = createOrOverwriteDataset<int>( gid, splitsName.c_str(), 2, dimsSplits,
      &(splits[0]) );
    if( err<0 ){
      cerr << "Dataset " << splitsName << " could not be written" << endl;
      return 1;
    }
  } else {
    if( H5LTfind_dataset(gid, splitsName.c_str()) ){
      herr_t status = H5Ldelete(gid, splitsName.c_str(), H5P_DEFAULT);
      if( status<0 ){
        cerr << "Error while unlinking dataset " << splitsName << endl;
        return 1;
      }
    }
  }
  const string deathsName = "Disappearances";
  if( !deaths.empty() ){
    hsize_t dimDeaths = deaths.size();
    herr_t err = createOrOverwriteDataset<int>( gid, deathsName.c_str(), 1,
      &dimDeaths, &(deaths[0]) );
    if( err<0 ){
      cerr << "Dataset " << deathsName << " could not be written" << endl;
      return 1;
    }
  } else {
    if( H5LTfind_dataset(gid, deathsName.c_str()) ){
      herr_t status = H5Ldelete( gid, deathsName.c_str(), H5P_DEFAULT );
      if( status<0 ){
        cerr << "Error while unlinking dataset " << deathsName << endl;
        return 1;
      }
    }
  }
  const string appearsName = "Appearances";
  if( !appears.empty() ){
    hsize_t dimAppears = appears.size();
   herr_t err = createOrOverwriteDataset<int>(gid, appearsName.c_str(), 1, 
      &dimAppears, &(appears[0]) );
    if ( err<0 ){
      cerr << "Dataset " << appearsName << " could not be written" << endl;
      return 1;
    }
  } else {
    if( H5LTfind_dataset( gid, appearsName.c_str() ) ){
      herr_t status = H5Ldelete( gid, appearsName.c_str(), H5P_DEFAULT );
      if( status<0 ){
        cerr << "Error while unlinking dataset " << appearsName << endl;
        return 1;
      }
    }
  }
  const string versionName = "Version";
  hsize_t dimScalar = 1;
  herr_t err = createOrOverwriteScalar<int>(gid, versionName.c_str(), &VERSION);
  if( err < 0 ){
    cerr << "Version tag could not be written" << endl;
    return 1;
  }
  const string subversionName = "Subversion";
  err = createOrOverwriteScalar<int>(gid, subversionName.c_str(),  &SUBVERSION);
  if( err < 0 ){
    cerr << "Subversion tag could not be written" << endl;
    return 1;
  } 
  const string solveFlagName = "Output code of solver";
  err = createOrOverwriteScalar<int>(gid, solveFlagName.c_str(), &solveFlag );
  if( err<0 ){
    cerr << solveFlagName << " could not be written" << endl;
    return 1;
  }
  const string solveStringName = "Explanation of solver output code";
  err = createOrOverwriteScalar<char>(gid, solveStringName.c_str(), solveString.c_str() );
  if( err<0 ){
    cerr << solveStringName << " could not be written" << endl;
    return 1;
  }
  const string totCostName = "Total costs";
  err = createOrOverwriteScalar<double>(gid, totCostName.c_str(), &totCost );
  if( err<0 ){
    cerr << "Dataset " << totCostName << " could not be written" << endl;
    return 1;
  }
  // Write the matching parameters
  const string knnName = "Number of nearest neighbors as possible children";
  herr_t err1 = createOrOverwriteScalar<int>(gid, knnName.c_str(), &pars.knn);
  const string maxDistName = "Maximum parent-child distance";
  herr_t err2 = createOrOverwriteScalar<double>(gid, maxDistName.c_str(), &pars.maxDist);
  const string cSplitName = "Cost of cell split";
  herr_t err3 = createOrOverwriteScalar<double>(gid, cSplitName.c_str(), &pars.cSplit);
  const string cDeathName = "Cost of cell death";
  herr_t err4 = createOrOverwriteScalar<double>(gid, cDeathName.c_str(), &pars.cDeath);
  const string cAppearName = "Cost of cell appearance";
  herr_t err5 = createOrOverwriteScalar<double>(gid, cAppearName.c_str(), &pars.cAppear);
  const string cFlagsName = "Feature extraction flags";
  herr_t err6 = createOrOverwriteScalar<unsigned long>(gid, cFlagsName.c_str(), 
    &pars.flagFeatures);
  const string xyzRatioName = "Ratio between lateral (XY) and axial (Z) resolution";
  herr_t err7 = createOrOverwriteScalar<double>(gid, xyzRatioName.c_str(), &pars.xyzRatio);
  const string solverNameName = "Name of MIP solver used for the association";
  herr_t err8 = createOrOverwriteScalar<char>(gid, solverNameName.c_str(), solverName.c_str());
  const string matchingTimeName = "CPU time for matching [sec]";
  herr_t err9 = createOrOverwriteScalar<double>(gid, matchingTimeName.c_str(), &matchingTime);
  bool errorOutput= (err1<0) || (err2<0) || (err3<0) || (err4<0) || (err5<0) ||
                    (err6<0) || (err7<0) || (err8<0) || (err9<0);
  if( errorOutput ){
    cerr << "Error while writing matching parameters" << endl;
  }
  err1 = H5Gclose( gid);
  err2 = H5Fclose( fid );
  if( (err1<0) || (err2<0) ){
    cerr << "Error while closing file or group" << endl;
  }
  return 0;
}

// Run a job on two subsequent files
// The output gives information about the successful completion
int runJob( const string& file1, const string& file2, const MatchingPars& pars,
            IlpSolver& solver ){
  int jobOutput = SUCCESSFUL;
  int nCells1 = 0, nCells2 = 0;
  int nFeats1 = 0, nFeats2 = 0;
  shared_array<double> feats1, feats2;
  FileReaderHci reader1( file1 );
  reader1.setFlag( pars.flagFeatures );
  int output = reader1.readInputs( nCells1, nFeats1, feats1 );
  if( output ){
    cerr << "Error reading coordinates from file " << file1 << endl;
    return READING_ERROR;
  }
  vector<int> validLabels1;
  reader1.getValidLabels( validLabels1 );
  FileReaderHci reader2( file2 );
  reader2.setFlag( pars.flagFeatures );
  output = reader2.readInputs( nCells2, nFeats2, feats2 );
  if( output ){
    cerr << "Error reading coordinates from file " << file2 << endl;
    return READING_ERROR;
  }
  vector<int> validLabels2;
  reader2.getValidLabels( validLabels2 );
  vector<int> moves, deaths, splits, appears;
  double totCost = 0;
  int solveFlag = 0;
  cout << "Computation of matches between " << file1 << " and " 
       << file2 << " started" << endl;
  string problemString = file1 + " -- " + file2;
  string solveString;
  clock_t startTime = clock();
  int outputMatch = computeMatch( nCells1, feats1, validLabels1, nCells2, feats2, 
                                  validLabels2, nFeats1, problemString, solver,
                                  moves, splits, deaths, appears, totCost,
                                  solveFlag, solveString, pars );
  clock_t endTime = clock();
  double timeForJob = ((double)(endTime - startTime)) / CLOCKS_PER_SEC;
  if( outputMatch ){
    cerr << "Error during main computation" << endl;
    jobOutput = COMPUTATION_ERROR;
  }
  output = writeOutputs( file2.c_str(), pars, moves, splits, deaths, appears,
                         totCost, solveFlag, solveString, timeForJob );
  if( output ){
    cerr << "Error while writing results to " << file2 << endl;
    jobOutput = WRITING_ERROR;
  }
  return jobOutput;
}

// Run the program
int main(int argc, char *argv[]){
  if( argc<2 ){
    printHelp();
    return 1;
  }
  vector<string> fileNames;
  if( fileExists( argv[1] ) ){
    try{
      std::ifstream input( argv[1] );
      std::copy( std::istream_iterator<string>(input), std::istream_iterator<string>(),
        std::back_inserter(fileNames));
    } catch( ... ){
      cerr << "Error reading data file names";
      return 1;
    }
  } else {
    cerr << "File with file names " << argv[1] << " does not exist: error";
    return 1;
  }
  size_t nFiles = fileNames.size();
  if( nFiles <= 1 ){
    cerr << "At least two files must be specified" << endl;
    return 1;
  }
  MatchingPars pars;
  int output = pars.readOptions( argc - 2, argv + 2 );
  if( output ){
    cerr << "Options wrongly specified";
    printHelp();
    return 1;
  }
  // initialize the MIP solver
#ifdef USE_CPLEX
  CplexSolver solver;
#else
  LpSolveSolver solver;
#endif
  for( size_t iF = 0; iF < nFiles-1; ++iF ){
    int output = runJob( fileNames[iF], fileNames[iF+1], pars, solver );
    cout << "Job no. " << iF+1 << " ended with output code " << output << endl;
  }
  cout << "Computation ended successfully" << endl;
  return 0;
}


