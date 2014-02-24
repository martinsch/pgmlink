# This module finds an installed Vigra package.
#
# It sets the following variables:
#  FFTW3_FOUND              - Set to false, or undefined, if fftw3 isn't found.
#  FFTW3_INCLUDE_DIR        - Fftw3 include directory.
#  FFTW3_LIBRARIES          - Fftw3 library files
FIND_PATH(FFTW3_INCLUDE_DIR fftw3.h PATHS /usr/include /usr/local/include ${CMAKE_INCLUDE_PATH} ${CMAKE_PREFIX_PATH}/include $ENV{FFTW3_ROOT}/include ENV CPLUS_INCLUDE_PATH)
FIND_LIBRARY(FFTW3_LIBRARIES fftw3 PATHS /usr/lib $ENV{FFTW3_ROOT}/lib ENV LD_LIBRARY_PATH ENV LIBRARY_PATH)

GET_FILENAME_COMPONENT(FFTW3_LIBRARY_PATH ${FFTW3_LIBRARIES} PATH)
SET( FFTW3_LIBRARY_DIR ${FFTW3_LIBRARY_PATH} CACHE PATH "Path to fftw3 library.")

# handle the QUIETLY and REQUIRED arguments and set FFTW3_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FFTW3 DEFAULT_MSG FFTW3_LIBRARIES FFTW3_INCLUDE_DIR)
IF(FFTW3_FOUND)
    IF (NOT Fftw3_FIND_QUIETLY)
      MESSAGE(STATUS "  > fftw3 includes: ${FFTW3_INCLUDE_DIR}")
      MESSAGE(STATUS "  > fftw3 libraries: ${FFTW3_LIBRARIES}")
    ENDIF()
ENDIF()

MARK_AS_ADVANCED( FFTW3_INCLUDE_DIR FFTW3_LIBRARIES)
