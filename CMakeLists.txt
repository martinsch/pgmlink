# our functions and macros
    # logging_level_to_define()
    # Gives the define corresponding to a global logging level.
    #
    # For example: The logging level 'INFO' corresponds to the define 'logINFO'.
    #
    # Parameters:
    #   LOGGING_LEVEL   STRING  One of the logging levels:
    #                           NO_LOGGING, ERROR, WARNING, INFO, DEBUG, DEBUG1, (...), DEBUG4
    #                           input parameter
    #
    #   DEFINE          STRING  output parameter
    # 
    MACRO(LOGGING_LEVEL_TO_DEFINE LOGGING_LEVEL DEFINE)
        IF(${LOGGING_LEVEL} STREQUAL "NO_LOGGING")
            SET(${DEFINE} "pgmlink::logNO_LOGGING")
        ELSEIF(${LOGGING_LEVEL} STREQUAL "ERROR")
            SET(${DEFINE} "pgmlink::logERROR")
        ELSEIF(${LOGGING_LEVEL} STREQUAL "WARNING")
            SET(${DEFINE} "pgmlink::logWARNING")
        ELSEIF(${LOGGING_LEVEL} STREQUAL "INFO")
            SET(${DEFINE} "pgmlink::logINFO")
        ELSEIF(${LOGGING_LEVEL} STREQUAL "DEBUG")
            SET(${DEFINE} "pgmlink::logDEBUG")
        ELSEIF(${LOGGING_LEVEL} STREQUAL "DEBUG1")
            SET(${DEFINE} "pgmlink::logDEBUG1")
        ELSEIF(${LOGGING_LEVEL} STREQUAL "DEBUG2")
            SET(${DEFINE} "pgmlink::logDEBUG2")
        ELSEIF(${LOGGING_LEVEL} STREQUAL "DEBUG3")
            SET(${DEFINE} "pgmlink::logDEBUG3")
        ELSEIF(${LOGGING_LEVEL} STREQUAL "DEBUG4")
            SET(${DEFINE} "pgmlink::logDEBUG4")
        ELSE(${LOGGING_LEVEL} STREQUAL "NO_LOGGING")
            MESSAGE(SEND_ERROR "Unknown LOGGING_LEVEL: ${LOGGING_LEVEL}. Default to INFO.")
            SET(${DEFINE} "pgmlink::logINFO")
        ENDIF(${LOGGING_LEVEL} STREQUAL "NO_LOGGING")
    ENDMACRO(LOGGING_LEVEL_TO_DEFINE)

###########################################################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(pgmlink)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_extensions/)
message( "\nConfiguring pgmlink:" )

# dependencies
find_package( ANN REQUIRED )
find_package( Cplex REQUIRED )
find_package( VIGRA REQUIRED )
find_package( Lemon REQUIRED )
find_package( Boost REQUIRED COMPONENTS serialization )
# find hdf5 for OSX
find_package( HDF5 REQUIRED COMPONENTS HL)
message(STATUS "  found: ${Boost_LIBRARIES}")

include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/include/ ${ANN_INCLUDE_DIR} ${CPLEX_INCLUDE_DIRS} ${VIGRA_INCLUDE_DIR} ${LEMON_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})

# CPLEX switch to be compatible with STL
ADD_DEFINITIONS(-DIL_STD)

## vigra defines
## assume that HDF5 is installed (for vigra random forest)
ADD_DEFINITIONS(-DHasHDF5)

# print out the dependencies
message(STATUS "Dependencies include dirs (you should check if they match the found libs):")
message(STATUS "  ANN: ${ANN_INCLUDE_DIR}")
message(STATUS "  Cplex: ${CPLEX_INCLUDE_DIR}")
message(STATUS "  VIGRA: ${VIGRA_INCLUDE_DIR}")
message(STATUS "  Lemon: ${LEMON_INCLUDE_DIR}")
message(STATUS "  Boost: ${Boost_INCLUDE_DIRS}")
message(STATUS "  HDF5:  ${HDF5_INCLUDE_DIRS}")

# hide these variables to have a cleaner interface
mark_as_advanced( Boost_DIR )

# build options

# logging level
set(LOGGING_LEVEL "INFO" CACHE STRING "Choose a global logging level: NO_LOGGING, ERROR, WARNING, INFO, DEBUG, DEBUG1, ..., DEBUG4")
# unit tests
set(WITH_TESTS "False" CACHE BOOL "Build tests.")
# build type and compiler options
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: Release Debug."
    FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if(NOT MSVC)
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wfatal-errors -std=c++98 -pedantic -pthread -fopenmp")
  set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Werror")
endif()

set(WITH_CHECKED_STL "True" CACHE BOOL "Build with checked stl for more robust code (may degrade performance; gcc only).")
if(WITH_CHECKED_STL)
  add_definitions(-D_GLIBCXX_DEBUG)
endif()

# global logging level
logging_level_to_define(LOGGING_LEVEL LOG_DEFINE)
add_definitions(-D FILELOG_MAX_LEVEL=${LOG_DEFINE})

# libpgmlink
file(GLOB LIB_SOURCES src/*.cpp)
add_library(pgmlink SHARED ${LIB_SOURCES})
# link with python for OSX
target_link_libraries(pgmlink ${ANN_LIBRARIES} ${Boost_LIBRARIES}  ${CPLEX_LIBRARIES} ${VIGRA_IMPEX_LIBRARY} ${LEMON_LIBRARIES} ${HDF5_LIBRARIES} mlpack fftw3 fftw3f)

# Install target pgmlink
install(TARGETS pgmlink
        DESTINATION lib)

# API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
else(DOXYGEN_FOUND)
  MESSAGE(WARNING "Doxygen not found. Won't add 'doc' build target.")
endif(DOXYGEN_FOUND)

# Tests
if(WITH_TESTS)
  enable_testing()
  add_subdirectory(tests/)
endif()

# Python
SET(WITH_PYTHON true CACHE BOOL "Build with python wrapper.")
if(WITH_PYTHON)
  add_subdirectory(python/)
endif()
