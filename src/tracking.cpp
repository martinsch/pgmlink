#include "common.h"
#include "tracking.h"
#include <ANN/ANN.h>
#include "file_reader.h"
#include <ctime>
#include "hdf5.h"
#include "hdf5_hl.h"
#include <boost/shared_ptr.hpp>

#define isnan( x ) ( (x) != (x) )

using std::vector;
using std::string;
using boost::shared_ptr;
using boost::shared_array;
using std::cerr;
using std::cout;
using std::endl;
 
namespace Tracking {
const size_t chunkSize = 50; // For creating a chunked data set

// the following size constant must be large enough to hold all strings we may want to write to an HDF5 file
const size_t stringBufferSize = 80;

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
    output = 1;
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
    solver.solve(nVars, nConstr, nNonZero, costs, rhs, matbeg, matcnt,
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


ReturnCode runJob( const string& file1, const string& file2, 
                   const MatchingPars& pars, IlpSolver& solver ){
  ReturnCode jobOutput = SUCCESSFUL;
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
  output = writeOutputs( file2, pars, moves, splits, deaths, appears,
                         totCost, solveFlag, solveString, solver.getName(), timeForJob );
  if( output ){
    cerr << "Error while writing results to " << file2 << endl;
    jobOutput = WRITING_ERROR;
  }
  return jobOutput;
}

int writeOutputs( const string& fileName, const MatchingPars& pars,
                  const vector<int>& moves, const vector<int>& splits, 
                  const vector<int>& deaths, const vector<int>& appears,
                  double totCost, int solveFlag, const string& solveString,
                  const string& solverName, double matchingTime ){
  hid_t fid = H5Fopen( fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
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
  const string revisionName = "Global SVN revision";
  err = createOrOverwriteScalar<char>(gid, revisionName.c_str(),REVISION_STRING.c_str());
  if( err < 0 ){
    cerr << revisionName <<  " could not be written" << endl;
  }
  const string dateName = "Date of last commit";
  err = createOrOverwriteScalar<char>(gid, dateName.c_str(), DATE_STRING.c_str());
  if( err<0 ){
    cerr << dateName << " could not be written" << endl;
  }
  const string authorName = "Author of last commit" ;
  err = createOrOverwriteScalar<char>(gid, authorName.c_str(), AUTHOR_STRING.c_str());
  if( err<0 ){
    cerr << authorName << " could not be written" << endl;
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

}
