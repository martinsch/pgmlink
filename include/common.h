#ifndef COMMON_H
#define COMMON_H

#include <string>
#include <sstream>
/* 
   Convenience functions that are needed in different parts of the tracking
   program
*/

namespace Tracking {

// Check whether a file exists
int fileExists( const std::string& fileName );

std::string inttostr(int number);


}

#endif /* COMMON_H */
