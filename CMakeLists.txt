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

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.6)
PROJECT(pgmlink)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_extensions/)
message( "\nConfiguring pgmlink:" )

# dependencies
find_package( ANN REQUIRED )
find_package( Cplex )
find_package( GUROBI )
find_package( VIGRA REQUIRED )
find_package( Lemon REQUIRED )
find_package( Boost REQUIRED COMPONENTS serialization program_options )
message(STATUS "  found: ${Boost_LIBRARIES}")
find_package( Armadillo REQUIRED )
find_package( Mlpack REQUIRED )
find_package( Opengm REQUIRED )
find_package( Xml2 REQUIRED )

# check which optimization package is found:
set(FORCE_GUROBI false CACHE BOOL "If GUROBI and CPLEX are found, by default we build against CPLEX. Setting this to true forces to build against GUROBI")
if(CPLEX_FOUND)
    set(OPTIMIZER_INCLUDE_DIRS ${CPLEX_INCLUDE_DIRS})
    set(OPTIMIZER_LIBRARIES ${CPLEX_LIBRARIES})

    if(GUROBI_FOUND AND FORCE_GUROBI)
        set(OPTIMIZER_INCLUDE_DIRS ${GUROBI_INCLUDE_DIR})
        set(OPTIMIZER_LIBRARIES ${GUROBI_LIBRARY} ${GUROBI_CXX_LIBRARY})
        add_definitions(-DWITH_GUROBI)
    endif()
else()
    if(GUROBI_FOUND)
        set(OPTIMIZER_INCLUDE_DIRS ${GUROBI_INCLUDE_DIR})
        set(OPTIMIZER_LIBRARIES ${GUROBI_LIBRARY} ${GUROBI_CXX_LIBRARY})
        add_definitions(-DWITH_GUROBI)
    else()
        message(ERROR "No optimizer found at all!")
    endif()
endif()

# hdf5
if(WIN32)
  # FindHDF5 is broken on Windows
  # it wrongly discovers the dlls instead of the lib files
  # see: www.cmake.org/Bug/bug_relationship_graph.php?bug_id=14111
  # therefore we set the variable manually
  find_library( HDF5_IMPORT_LIB NAMES hdf5 hdf5dll )
  find_library( HDF5_HL_IMPORT_LIB NAMES hdf5_hl hdf5_hldll )
  set( HDF5_LIBRARIES ${HDF5_IMPORT_LIB} ${HDF5_HL_IMPORT_LIB} )
  message(STATUS "HDF5: autodiscovery is broken on WIN32; using global search paths")
  message(STATUS "HDF5: found ${HDF5_LIBRARIES}")
else()
  find_package( HDF5 REQUIRED COMPONENTS C HL)
endif()


include_directories(${PROJECT_SOURCE_DIR}/include/)
# include external headers as system includes so we do not have to cope with their warnings
include_directories(SYSTEM ${ANN_INCLUDE_DIR} ${OPTIMIZER_INCLUDE_DIRS} ${VIGRA_INCLUDE_DIR} ${LEMON_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${Opengm_INCLUDE_DIR} ${Xml2_INCLUDE_DIR})

# CPLEX switch to be compatible with STL
ADD_DEFINITIONS(-DIL_STD)

if(WIN32)
  add_definitions(-DBOOST_ALL_DYN_LINK)
  if("${Boost_LIBRARIES}" MATCHES "boost_serialization\\.lib")
    ADD_DEFINITIONS(-DBOOST_AUTO_LINK_NOMANGLE)
  endif()
endif()

## vigra defines
## assume that HDF5 is installed (for vigra random forest)
ADD_DEFINITIONS(-DHasHDF5)

# opengm defines: do not include LPDef symbols everywhere
ADD_DEFINITIONS(-DOPENGM_LPDEF_NO_SYMBOLS)

# print out the dependencies
message(STATUS "Dependencies include dirs (you should check if they match the found libs):")
message(STATUS "  ANN: ${ANN_INCLUDE_DIR}")
message(STATUS "  Optimizer: ${OPTIMIZER_INCLUDE_DIR}")
message(STATUS "  VIGRA: ${VIGRA_INCLUDE_DIR}")
message(STATUS "  Lemon: ${LEMON_INCLUDE_DIR}")
message(STATUS "  Boost: ${Boost_INCLUDE_DIRS}")
message(STATUS "  HDF5:  ${HDF5_INCLUDE_DIRS}")
message(STATUS "  Opengm: ${Opengm_INCLUDE_DIR}")
message(STATUS "  libxml2:  ${Xml2_INCLUDE_DIR}")

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

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -std=c++98 -pedantic -pthread -fopenmp ${CMAKE_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 ${CMAKE_CXX_FLAGS_DEBUG}")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  add_definitions("-stdlib=libstdc++")
  set(CMAKE_CXX_FLAGS "-std=c++98 -ftemplate-depth=512 ${CMAKE_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 ${CMAKE_CXX_FLAGS_DEBUG}")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

IF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(WITH_CHECKED_STL "True" CACHE BOOL "Build with checked stl for more robust code (may degrade performance; gcc only).")
  if(WITH_CHECKED_STL)
    add_definitions(-D_GLIBCXX_DEBUG)
  endif()
ENDIF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 

# global logging level
logging_level_to_define(LOGGING_LEVEL LOG_DEFINE)
add_definitions(-D FILELOG_MAX_LEVEL=${LOG_DEFINE})

# libpgmlink
include( GenerateExportHeader )
## only activate symbol export on Windows
if(NOT WIN32)
  set(USE_COMPILER_HIDDEN_VISIBILITY OFF)
endif()

if(WIN32)
  add_compiler_export_flags() # enable dynamic symbol export/import for shared libs (on Windows)
endif()

file(GLOB LIB_SOURCES src/*.cpp src/feature_calculator/*.cpp)
file(GLOB_RECURSE HEADERS include/*.h*)

add_library(pgmlink SHARED ${LIB_SOURCES} ${HEADERS})



file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/pgmlink )
if(WIN32)
  generate_export_header( pgmlink ) # generates pgmlink_export.h
  file( RENAME ${CMAKE_CURRENT_BINARY_DIR}/pgmlink_export.h  ${CMAKE_CURRENT_BINARY_DIR}/include/pgmlink/pgmlink_export.h )
else(WIN32)
  configure_file(${CMAKE_SOURCE_DIR}/include/pgmlink/pgmlink_export.h.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/include/pgmlink/pgmlink_export.h
    COPYONLY)
endif(WIN32)


include_directories( ${CMAKE_CURRENT_BINARY_DIR}/include )

target_link_libraries(pgmlink ${ANN_LIBRARIES} ${Boost_LIBRARIES}  ${OPTIMIZER_LIBRARIES} ${VIGRA_IMPEX_LIBRARY} ${LEMON_LIBRARIES} ${HDF5_LIBRARIES} ${Mlpack_LIBRARIES} ${Armadillo_LIBRARIES})

# Install target pgmlink
install(TARGETS pgmlink
	ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)

# install headers
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/"
        DESTINATION include
        PATTERN ".h"
        PATTERN ".hxx"
        PATTERN ".git" EXCLUDE)
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include/"
        DESTINATION include
        PATTERN ".h"
        PATTERN ".git" EXCLUDE)

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
  find_package( Vigranumpy REQUIRED PATHS ${Vigranumpy_DIR} )
  SET( Vigranumpy_DIR ${Vigranumpy_DIR}
    CACHE PATH "Directory containing VigranumpyConfig.cmake"
    FORCE )
  include_directories(SYSTEM ${VIGRANUMPY_INCLUDE_DIRS})
  add_subdirectory(python/)
endif()
