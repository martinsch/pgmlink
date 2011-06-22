#include "common.h"
#include <sys/stat.h>

namespace Tracking {

// Check whether the file exists
int fileExists( const std::string& fileName ) {
  int output=0;
  if( !fileName.empty() ){
    struct stat buf;
    output = (stat( fileName.c_str(), &buf ) == 0);
  }
  return output;
}


std::string inttostr(int number)
{
    std::stringstream ss; //create a stringstream
    ss << number;         //add number to the stream
    return ss.str();      //return a string with the contents of the stream
}

}

