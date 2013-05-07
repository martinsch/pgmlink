# This module finds an installed Lemon package.
#
# It sets the following variables:
#  LEMON_FOUND              - Set to false, or undefined, if lemon isn't found.
#  LEMON_INCLUDE_DIR        - Lemon include directory.
#  LEMON_LIBRARIES          - Lemon library files
FIND_PATH(LEMON_INCLUDE_DIR lemon/config.h PATHS /usr/include /usr/local/include ${CMAKE_INCLUDE_PATH} ${CMAKE_PREFIX_PATH}/include $ENV{LEMON_ROOT}/include ENV CPLUS_INCLUDE_PATH)
FIND_LIBRARY(LEMON_LIBRARIES 
  NAMES emon lemon 
  PATHS $ENV{LEMON_ROOT}/src/impex $ENV{LEMON_ROOT}/lib ENV LD_LIBRARY_PATH ENV LIBRARY_PATH
)

GET_FILENAME_COMPONENT(LEMON_LIBRARY_PATH ${LEMON_LIBRARIES} PATH)
SET( LEMON_LIBRARY_DIR ${LEMON_LIBRARY_PATH} CACHE PATH "Path to lemon library.")

# handle the QUIETLY and REQUIRED arguments and set LEMON_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LEMON DEFAULT_MSG LEMON_LIBRARIES LEMON_INCLUDE_DIR)

MARK_AS_ADVANCED( LEMON_INCLUDE_DIR LEMON_LIBRARIES LEMON_LIBRARY_DIR )
