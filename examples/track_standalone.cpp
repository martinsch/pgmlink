#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <boost/shared_ptr.hpp>
#include "tracking.h"
#include "common.h"
#ifdef USE_CPLEX
#include "cplex_solver.h"
#else
#include "lp_solve_solver.h"
#endif

using namespace Tracking;
using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::ostringstream;
using boost::shared_array;
using boost::shared_ptr;
 
// Print the usage instructions for this program 
void printHelp() {
  cout << "track_standalone : Simple tracking of presegmented cells\n"
       << "Mercurial changeset: " << REVISION_STRING << "\n"
       << "Version " << VERSION << "." << SUBVERSION  << "\n"
       << "Last change: " << DATE_STRING << "\n"
       << "Author of last change: " << AUTHOR_STRING << "\n"
       << "Usage: track_standalone TASKFILE [OPTIONS]\n"
       << "TASKFILE : text file with the names of all files on which the program"
          " shall operate\n"
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


