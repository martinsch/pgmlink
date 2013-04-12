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


##
# dependencies
##
    find_package( VIGRA REQUIRED )
    find_package( Lemon REQUIRED )
    find_package( HDF5 REQUIRED )
    find_package( FFTW3 REQUIRED )
#    find_package( PythonLibs REQUIRED )
#    find_package(Vigranumpy REQUIRED PATHS ${Vigranumpy_DIR})

    include_directories(${FFTW3_INCLUDE_DIR} ${VIGRA_INCLUDE_DIR} ${LEMON_INCLUDE_DIR} ${HDF5_INCLUDE_DIRS} /usr/include/python2.7/ )
    link_directories( ${FFTW3_LIBRARY_DIR} ${VIGRA_IMPEX_LIBRARY_DIR} ${LEMON_LIBRARY_DIR} ${HDF5_LIBRARY_DIRS} )
##
# build options
##
#    SET(Vigranumpy_DIR ${Vigranumpy_DIR}
#      CACHE PATH "Directory containing VigranumpyConfig.cmake" FORCE )
#    message(STATUS "Vigranumpy found at: ${Vigranumpy_CONFIG}") 

    # logging level
    SET(LOGGING_LEVEL "INFO" CACHE STRING "Choose a global logging level: NO_LOGGING, ERROR, WARNING, INFO, DEBUG, DEBUG1, ..., DEBUG4")

    # unit tests
    SET(WITH_TESTS "False" CACHE BOOL "Build tests.")

##
# build type and compiler options
##
    IF(NOT CMAKE_BUILD_TYPE)
      SET(CMAKE_BUILD_TYPE Release CACHE STRING
	"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
	FORCE)
    ENDIF(NOT CMAKE_BUILD_TYPE)

    SET(CMAKE_CXX_FLAGS "-Wall -Wextra -Wfatal-errors -std=c++98 -pedantic -pthread -fopenmp")
    SET(CMAKE_CXX_FLAGS_DEBUG "-Wno-long-long -g -O0")
    SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -Werror")

    SET(WITH_CHECKED_STL "True" CACHE BOOL "Build with checked libg++ for more robust code (may degrade performance).")
    IF(WITH_CHECKED_STL)
      ADD_DEFINITIONS(-D_GLIBCXX_DEBUG)
    ENDIF()
    
    #CPLEX switch to be compatible with STL
    ADD_DEFINITIONS(-DIL_STD)

    # assume by default that HDF5 is installed. Set definition for vigra random forest
    ADD_DEFINITIONS(-DHasHDF5)

##
# global logging level
##
    LOGGING_LEVEL_TO_DEFINE(LOGGING_LEVEL LOG_DEFINE)
    ADD_DEFINITIONS(-D FILELOG_MAX_LEVEL=${LOG_DEFINE})

INCLUDE_DIRECTORIES(
  "${PROJECT_SOURCE_DIR}/include/" 
)

FILE(GLOB LIB_SOURCES src/*.cpp)

find_package( Boost REQUIRED COMPONENTS serialization )

ADD_LIBRARY(pgmlink SHARED ${LIB_SOURCES})
#TARGET_LINK_LIBRARIES(pgmlink ${LIBS} ${VIGRA_NUMPY_CORE_LIBRARY} vigraimpex boost_python ${PYTHON_LIBRARIES} ${FFTW3_LIBRARIES} fftw3 fftw3f fftw3_threads ${HD5_LIBRARIES} ${HDF5_HL_LIBRARIES} ann cplex ilocplex pthread concert ${VIGRA_IMPEX_LIBRARY} ${LEMON_LIBRARIES} )
TARGET_LINK_LIBRARIES(pgmlink ${LIBS} boost_python ${PYTHON_LIBRARIES} python2.7 ${FFTW3_LIBRARIES} fftw3 fftw3f fftw3_threads ${HD5_LIBRARIES} ${HDF5_HL_LIBRARIES} ann cplex ilocplex pthread concert ${VIGRA_IMPEX_LIBRARY} ${LEMON_LIBRARIES} mlpack)


### Always use full RPATH

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

##
# Tests
if(WITH_TESTS)
  enable_testing()
  add_subdirectory(tests/)
endif()

##
# Python
SET(WITH_PYTHON true CACHE BOOL "build with python wrapper")
if(WITH_PYTHON)
  add_subdirectory(python/)
endif()

##
# API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)
add_custom_target(doc
${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)

##
# Install target pgmlink
INSTALL(TARGETS pgmlink
        DESTINATION lib)

