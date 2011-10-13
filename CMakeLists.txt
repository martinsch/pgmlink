##
# our functions and macros
##
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
            SET(${DEFINE} "Tracking::logNO_LOGGING")
        ELSEIF(${LOGGING_LEVEL} STREQUAL "ERROR")
            SET(${DEFINE} "Tracking::logERROR")
        ELSEIF(${LOGGING_LEVEL} STREQUAL "WARNING")
            SET(${DEFINE} "Tracking::logWARNING")
        ELSEIF(${LOGGING_LEVEL} STREQUAL "INFO")
            SET(${DEFINE} "Tracking::logINFO")
        ELSEIF(${LOGGING_LEVEL} STREQUAL "DEBUG")
            SET(${DEFINE} "Tracking::logDEBUG")
        ELSEIF(${LOGGING_LEVEL} STREQUAL "DEBUG1")
            SET(${DEFINE} "Tracking::logDEBUG1")
        ELSEIF(${LOGGING_LEVEL} STREQUAL "DEBUG2")
            SET(${DEFINE} "Tracking::logDEBUG2")
        ELSEIF(${LOGGING_LEVEL} STREQUAL "DEBUG3")
            SET(${DEFINE} "Tracking::logDEBUG3")
        ELSEIF(${LOGGING_LEVEL} STREQUAL "DEBUG4")
            SET(${DEFINE} "Tracking::logDEBUG4")
        ELSE(${LOGGING_LEVEL} STREQUAL "NO_LOGGING")
            MESSAGE(SEND_ERROR "Unknown LOGGING_LEVEL: ${LOGGING_LEVEL}. Default to INFO.")
            SET(${DEFINE} "Tracking::logINFO")
        ENDIF(${LOGGING_LEVEL} STREQUAL "NO_LOGGING")
    ENDMACRO(LOGGING_LEVEL_TO_DEFINE)

###########################################################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(tracking)

SET(VERSION 1)
SET(SUBVERSION 5)
SET(ANN_VERSION 1.1.2)
SET(LP_SOLVE_VERSION 5)
SET(LP_SOLVE_SUBVERSION 5)
SET(ARCHITECTURE ux64)

##
# build options
##
    # logging level
    SET(LOGGING_LEVEL "INFO" CACHE STRING "Choose a global logging level: NO_LOGGING, ERROR, WARNING, INFO, DEBUG, DEBUG1, ..., DEBUG4")

    # unit tests
    set(WITH_TESTS "False" CACHE BOOL "Build tests.")

##
# build type and compiler options
##
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -std=c++98 -pedantic")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Werror")

    set(WITH_CHECKED_STL "True" CACHE BOOL "Build with checked libg++ for more robust code (may degrade performance).")
    if(WITH_CHECKED_STL)
      ADD_DEFINITIONS(-D_GLIBCXX_DEBUG)
    endif()
    
    ADD_DEFINITIONS(-DUSE_CPLEX)
    #CPLEX switch to be compatible with STL
    ADD_DEFINITIONS(-DIL_STD)

    # assume by default that HDF5 is installed. Set definition for vigra random forest
    ADD_DEFINITIONS(-DHasHDF5)

##
# global logging level
##
    LOGGING_LEVEL_TO_DEFINE(LOGGING_LEVEL LOG_DEFINE)
    ADD_DEFINITIONS(-D FILELOG_MAX_LEVEL=${LOG_DEFINE})


# pass information about the revision number to the compilation
#FIND_PROGRAM(HG_EXECUTABLE hg DOC "mercurial command line client")
#EXECUTE_PROCESS(COMMAND ${HG_EXECUTABLE} --debug tip OUTPUT_VARIABLE hginfo OUTPUT_STRIP_TRAILING_WHITESPACE)
#STRING(REGEX REPLACE "^(.*\n)?changeset:([ ]+)([^\n]+).*" "\\3" ProjectRevision ${hginfo})
#STRING(REGEX REPLACE "^(.*\n)?user:([ ]+)([^\n]+).*" "\\3" LastChangedAuthor ${hginfo})
#STRING(REGEX REPLACE "^(.*\n)?date:([ ]+)([^\n]+).*" "\\3" LastChangedDate ${hginfo})

#CONFIGURE_FILE(
#  "${PROJECT_SOURCE_DIR}/include/tracking.h.in"
#  "${PROJECT_BINARY_DIR}/include/tracking.h"
#)

INCLUDE_DIRECTORIES(
  "${PROJECT_SOURCE_DIR}/include/" 
)
#LINK_DIRECTORIES( 
#  ${PROJECT_SOURCE_DIR}/../3rd_party/lib
#)




FILE(GLOB LIB_SOURCES src/*.cpp)

ADD_LIBRARY(track_lib SHARED ${LIB_SOURCES})
TARGET_LINK_LIBRARIES(track_lib hdf5 hdf5_hl ann cplex ilocplex pthread concert vigraimpex emon)
#TARGET_LINK_LIBRARIES(track_lib lpsolve${LP_SOLVE_VERSION}${LP_SOLVE_SUBVERSION})

#ADD_EXECUTABLE(track_standalone examples/track_standalone.cpp)
#TARGET_LINK_LIBRARIES(track_standalone track_lib vigraimpex concert)

# Tests
if(WITH_TESTS)
  enable_testing()
  add_subdirectory(tests/)
endif()

# Python
set(WITH_PYTHON true CACHE BOOL "build with python wrapper")
if(WITH_PYTHON)
  add_subdirectory(python/)
endif()

# Binaries
set(WITH_BINARIES true CACHE BOOL "binaries using the library")
if(WITH_BINARIES)
  add_subdirectory(bin/)
endif()
