#include "lp_lib.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "ANN.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include <sys/stat.h>

using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using boost::shared_array;
using boost::shared_ptr;

int fileExists( char* fileName ){
  struct stat buf;
  int i = stat( fileName, &buf );
  return i==0;
}

/* C/C++ implementation of the module track2.m (see comments in this file)

   Track cells in a greedy fashion (time slice by time slice) by 
   minimizing a simple energy function (mainly based on distances) via
   Integer Linear Programming.

   Implementation is geared towards trivial (data-parallel) parallelization.

   Requires the following 3rd party libraries:
   - LP_SOLVE for mixed-integer programming 
     http://lpsolve.sourceforge.net/5.5/
   - HDF5 for input and output data storage
     http://www.hdfgroup.org/HDF5/
   - ANN for efficient nearest-neighbor search
     http://www.cs.umd.edu/~mount/ANN/

   Frederik Orlando Kaster, 2010 <frederik.kaster@iwr.uni-heidelberg.de>
 */

// Parameters for the matching process 
class MatchingPars {
public:
  int knn;
  double maxDist;
  double cSplit;
  double cDeath;
  double cAppear;
  MatchingPars();
  int readOptions(int nOpts, char* opts[]);
};

// Define default parameters in the constructor
MatchingPars::MatchingPars() : knn(6),
                               maxDist(25),
                               cSplit(100),
                               cDeath(1000),
                               cAppear(10000) {
}

/* Read parameters from the input options: see function printHelp() for
   explanations
*/
int MatchingPars::readOptions(int nOpts, char* opts[]){
  if( nOpts % 2 ){
    cerr << "Options must be passed as key-value pairs" << endl;
    return 1;
  }
  nOpts/=2;
  const string allowedKeys("kmsda");
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
      }
      if( !valueReader.good() ){
        cerr << "Invalid value for option no. " << iOpt << endl;
        return 1;
      }
    }
  }
  return 0;
} 

/* Print the usage instructions for this program */
void printHelp() {
  cout << "TRACK2 : Simple tracking of presegmented cells\n"
       << "Version 1.0, 2010/02/01, Frederik Kaster\n\n"
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
       << endl; 
}

/* Read the input coordinates from a specified HDF5 file 

   Upon successful completion, nCells contains the number of cells and
   coords contains the 3*nCells center-of-mass coordinates (first all the
   X coordinates stacked behind each other, then all the Y coordinates, then
   all the Z coordinates).
*/

int readInputs( const char* fileName, int& nCells, 
                shared_array<double>& coords ){
  const string dsetName = "/COM Coordinates";
  hid_t fid = H5Fopen( fileName, H5F_ACC_RDONLY, H5P_DEFAULT );
  if( fid<0 ){
    cerr << "File " << fileName << " cannot be opened" << endl;
    return 1;
  }
  int dataRank = 0;
  herr_t err = H5LTget_dataset_ndims( fid, dsetName.c_str(), &dataRank );
  if( err<0 ){
    cerr << "Dataset '" + dsetName + "' does not exist in " << fileName << endl;
    return 1;
  }
  if( dataRank!=2 ){
    cerr << "Dataset '" + dsetName + "' in " << fileName << " has wrong rank"
         << endl;
    return 1;
  }
  hsize_t dataDims[2] = { 0, 0 };
  H5T_class_t dataType = H5T_FLOAT;
  size_t dataTypeSize = 0;
  err = H5LTget_dataset_info( fid, dsetName.c_str(), dataDims, &dataType,
                              &dataTypeSize );
  if( err<0 ){
    cerr << "Error while getting properties of '" + dsetName + "' in file "
         << fileName << endl;
    return 1;
  }
  if( dataDims[0]!=3 ){
    cerr << "Data set " + dsetName + " in file " + fileName + " of wrong "
      "size" << endl;
    return 1;
  }
  if( dataType!=H5T_FLOAT ){
    cerr << "Data set " + dsetName + " in file " + fileName + " of wrong type"
   << endl;
    return 1;
  }
  nCells = dataDims[1];
  coords.reset( new double[3*nCells] );
  err = H5LTread_dataset_double( fid, dsetName.c_str(), coords.get() );
  if( err<0 ){
    cerr << "Error reading data from "+dsetName+" in file "+fileName << endl;
    return 1;
  }
  err = H5Fclose( fid );
  if( err < 0 ){
    cerr << "Error while closing file" << endl;
    return 1;
  }
  return 0;
}

/* Compute the optimum matching by optimizing a linear integer program.

   The inputs are usually acquired by the readInputs() function and the
   MatchingPars::readOptions() function.
   Upon successful completion, nMoves / nSplits / nDeaths / nAppears
   denote the numbers of motion / mitosis / death or disappearance /
   appearance events.

   moves : nMoves * 2 matrix with cell indices in the first and second
           slice
   splits : nSplits * 3 matrix, with the cell indices of the mother cells in
            the first row, and the indices of the daughter cells in the 
            second and third row
   deaths : nDeaths * 1 matrix, with the indices of the disappearing cells
            from the first time slice
   appears : nAppears * 1 matrix, with the indices of the appearing cells 
             from the second time slice

   All matrices are vectorized in row-major order (as usual in C).
   solveFlag shows the output flag of the LP_SOLVE operation, see

   http://lpsolve.sourceforge.net/5.5/solve.htm

   for the interpretation of the various integral return values.

   totCost holds the total costs of the optimum matching.
*/

int computeMatch( int nCells1, shared_array<double> coords1,
                  int nCells2, shared_array<double> coords2,
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
    const int dim = 3; // dimensionality of search space
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
    pts2 = annAllocPts( nCells2, dim );
    if( pts2==NULL ){
      throw "Allocation of points for kd-tree failed";
    }
    for( int i=0; i<nCells2; ++i ){
      ANNpoint cp = pts2[i];
      for( int j=0; j<dim; ++j ){
        cp[j] = coords2[j*nCells2+i];
      }
    }
    currDists = new ANNdist[pars.knn];
    shared_ptr<ANNkd_tree> kdTree( new ANNkd_tree( pts2, nCells2, dim ) );
    currQ = annAllocPt(dim);
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
    // MOTHER_INDEX DAUGHTER_INDEX 0 (for a move)
    // MOTHER_INDEX DAUGHTER1_INDEX DAUGHTER2_INDEX (for a split)
    vector<int> variables;
    for( int i=0; i<nCells1; ++i ){
      for( int j=0; j<dim; ++j ){
        currQ[j]=coords1[j*nCells1+i];
      } 
      int totPts = kdTree->annkFRSearch( currQ, pars.maxDist*pars.maxDist,
	 pars.knn, currNn, currDists, eps );
      totPts = (totPts < pars.knn) ? totPts : pars.knn;
      // add a column for every potential move
      for( int k=0; k<totPts; ++k ){
        moveColData[0] = currDists[k]-pars.cDeath-pars.cAppear;
        moveColData[1] = -1.;
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
        variables.push_back( i+1 );
        variables.push_back( currNn[k]+1 );
        variables.push_back( 0 );
        if( !out ){
          throw "Problem while constraining variable to be binary";
        }
      }
      // add a constraint column for every potential split
      for( int k=0; k<totPts; ++k ){
        for( int l=k+1; l<totPts; ++l ){
          splitColData[0] = currDists[k]+currDists[l]+pars.cSplit
				-pars.cDeath-2*pars.cAppear;
          splitColData[1] = -1.;
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
          variables.push_back( i+1 );
          variables.push_back( currNn[k]+1 );
          variables.push_back( currNn[l]+1 );
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
      rhs[i] = -1;
      out = set_constr_type( lp, i, GE );
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
      if( variables[3*i+2] == 0 ){ // move
        moves.push_back(variables[3*i]);
        moves.push_back(variables[3*i+1]);
        accountCells1[variables[3*i]-1]=true;
        accountCells2[variables[3*i+1]-1]=true;
      } else { // split
        splits.push_back(variables[3*i]);
        splits.push_back(variables[3*i+1]);
        splits.push_back(variables[3*i+2]);
        accountCells1[variables[3*i]-1]=true;
        accountCells2[variables[3*i+1]-1]=true;
        accountCells2[variables[3*i+2]-1]=true;
      }
    }
    // the cells unaccounted for must be deaths or appearances
    for( int i=0; i<nCells1; ++i ){
      if( !accountCells1[i] ){
        deaths.push_back( i+1 );
      }
    }
    for( int i=0; i<nCells2; ++i ){
      if( !accountCells2[i] ){
        appears.push_back( i+1 );
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

int writeOutputs( const char* fileName, 
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
  const string solveFlagName = "/Output flag of LP_SOLVE";
  hsize_t dimScalar = 1;
  herr_t err = H5LTmake_dataset_int(fid, solveFlagName.c_str(), 1, &dimScalar, 
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
  if( pars.readOptions( argc-4, argv+4) ){
    cerr << "Options wrongly specified\n";
    printHelp();
    return 2;
  }
  int nCells1 = 0, nCells2 = 0;
  shared_array<double> coords1, coords2;
  if( readInputs( argv[1], nCells1, coords1 ) ){
    cerr << "Error reading coordinates from " << argv[1] << endl;
    return 5;
  }
  if( readInputs( argv[2], nCells2, coords2 ) ){
    cerr << "Error reading coordinates from " << argv[2] << endl;
    return 6;
  }
  vector<int> moves, deaths, splits, appears;
  double totCost = 0;
  int solveFlag = 0;
  cout << "Computation of matches between " << argv[1] << " and " 
       << argv[2] << " started" << endl;
  if( computeMatch( nCells1, coords1, nCells2, coords2, moves, splits,
                    deaths, appears, totCost, solveFlag, pars ) ){
    cerr << "Error during main computation" << endl;
    return 7;
  }
  if( writeOutputs( argv[3], moves, splits, deaths, appears,
                    totCost, solveFlag ) ){
    cerr << "Error while writing results to " << argv[3] << endl;
    return 8;
  }
  cout << "Computation ended successfully" << endl;
  return 0;
}

