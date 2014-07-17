# This module finds an installed armadilo package.
#
# It sets the following variables:
#  Xml2_FOUND              - Set to false, or undefined, if lemon isn't found.
#  Xml2_INCLUDE_DIR        - Lemon include directory.
#  Xml2_LIBRARIES          - Lemon library files
FIND_PATH(Xml2_INCLUDE_DIR libxml/xmlversion.h PATHS /usr/include/libxml2 /usr/local/include/libxml2 ${CMAKE_INCLUDE_PATH} ${CMAKE_PREFIX_PATH}/include/libxml2 $ENV{Xml2_ROOT}/include/libxml2 ENV CPLUS_INCLUDE_PATH)
FIND_LIBRARY(Xml2_LIBRARIES 
  NAMES xml2 libxml2
  PATHS $ENV{Xml2_ROOT}/lib ${CMAKE_PREFIX_PATH}/lib ENV LD_LIBRARY_PATH ENV LIBRARY_PATH
)

GET_FILENAME_COMPONENT(Xml2_LIBRARY_PATH ${Xml2_LIBRARIES} PATH)
SET( Xml2_LIBRARY_DIR ${Xml2_LIBRARY_PATH} CACHE PATH "Path to lemon library.")

# handle the QUIETLY and REQUIRED arguments and set Xml2_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Xml2 DEFAULT_MSG Xml2_LIBRARIES Xml2_INCLUDE_DIR)

MARK_AS_ADVANCED( Xml2_INCLUDE_DIR Xml2_LIBRARIES Xml2_LIBRARY_DIR )
