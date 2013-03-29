# This module finds cplex.
#
# It sets the following variables:
#  CPLEX_FOUND              - Set to false, or undefined, if lemon isn't found.
#  CPLEX_INCLUDE_DIR        - include directory.
#  CPLEX_LIBRARIES          - library files
FIND_PATH(CPLEX_INCLUDE_DIR ilcplex/cplex.h PATHS /usr/include /usr/local/include ${CMAKE_INCLUDE_PATH} ${CMAKE_PREFIX_PATH}/include ENV C_INCLUDE_PATH ENV CPLUS_INCLUDE_PATH)
FIND_LIBRARY(cplex_LIBRARY cplex PATHS ENV LD_LIBRARY_PATH ENV LIBRARY_PATH)
FIND_LIBRARY(cplex_ilocplex_LIBRARY ilocplex PATHS ENV LD_LIBRARY_PATH ENV LIBRARY_PATH)
FIND_LIBRARY(cplex_concert_LIBRARY concert PATHS ENV LD_LIBRARY_PATH ENV LIBRARY_PATH)

if(cplex_LIBRARY AND cplex_ilocplex_LIBRARY AND cplex_ilocplex_LIBRARY)
  set(CPLEX_LIBRARIES ${cplex_concert_LIBRARY} ${cplex_ilocplex_LIBRARY} ${cplex_LIBRARY} )
else(cplex_LIBRARY AND cplex_ilocplex_LIBRARY AND cplex_ilocplex_LIBRARY)
  set(CPLEX_LIBRARIES FALSE)
endif(cplex_LIBRARY AND cplex_ilocplex_LIBRARY AND cplex_ilocplex_LIBRARY)

GET_FILENAME_COMPONENT(CPLEX_LIBRARY_PATH $CPLEX_LIBRARIES} PATH)
SET( CPLEX_LIBRARY_DIR ${CPLEX_LIBRARY_PATH} CACHE PATH "Path to cplex libraries.")

# handle the QUIETLY and REQUIRED arguments and set CPLEX_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CPLEX DEFAULT_MSG CPLEX_LIBRARIES CPLEX_INCLUDE_DIR)

MARK_AS_ADVANCED( CPLEX_INCLUDE_DIR CPLEX_LIBRARIES CPLEX_LIBRARY_DIR cplex_LIBRARY cplex_ilocplex_LIBRARY cplex_concert_LIBRARY )
