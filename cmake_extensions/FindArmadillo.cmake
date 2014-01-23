# This module finds an installed armadilo package.
#
# It sets the following variables:
#  Armadillo_FOUND              - Set to false, or undefined, if lemon isn't found.
#  Armadillo_INCLUDE_DIR        - Lemon include directory.
#  Armadillo_LIBRARIES          - Lemon library files
FIND_PATH(Armadillo_INCLUDE_DIR armadillo PATHS /usr/include /usr/local/include ${CMAKE_INCLUDE_PATH} ${CMAKE_PREFIX_PATH}/include $ENV{Armadillo_ROOT}/include ENV CPLUS_INCLUDE_PATH)
FIND_LIBRARY(Armadillo_LIBRARIES 
  NAMES armadillo 
  PATHS $ENV{Armadillo_ROOT}/src/impex $ENV{Armadillo_ROOT}/lib ENV LD_LIBRARY_PATH ENV LIBRARY_PATH
)

GET_FILENAME_COMPONENT(Armadillo_LIBRARY_PATH ${Armadillo_LIBRARIES} PATH)
SET( Armadillo_LIBRARY_DIR ${Armadillo_LIBRARY_PATH} CACHE PATH "Path to lemon library.")

# handle the QUIETLY and REQUIRED arguments and set Armadillo_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Armadillo DEFAULT_MSG Armadillo_LIBRARIES Armadillo_INCLUDE_DIR)

MARK_AS_ADVANCED( Armadillo_INCLUDE_DIR Armadillo_LIBRARIES Armadillo_LIBRARY_DIR )
