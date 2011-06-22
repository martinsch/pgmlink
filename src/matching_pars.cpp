#include "matching_pars.h"
#include <string>
#include <sstream>
#include <iostream>
using std::string;
using std::cerr;
using std::endl;

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
