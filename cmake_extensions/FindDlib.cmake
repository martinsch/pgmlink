# This module finds an installed armadilo package.
#
# It sets the following variables:
#  Dlib_FOUND              - Set to false, or undefined, if lemon isn't found.
#  Dlib_INCLUDE_DIR        - Dlib include directory.
FIND_PATH(Dlib_INCLUDE_DIR dlib/config.h PATHS /usr/include /usr/local/include ${CMAKE_INCLUDE_PATH} ${CMAKE_PREFIX_PATH}/include $ENV{Dlib_ROOT}/include ENV CPLUS_INCLUDE_PATH)

# handle the QUIETLY and REQUIRED arguments and set Dlib_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Dlib DEFAULT_MSG Dlib_INCLUDE_DIR)

MARK_AS_ADVANCED( Dlib_INCLUDE_DIR )
