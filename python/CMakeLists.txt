
cmake_minimum_required(VERSION 2.8)

##
# dependencies
##
find_package( Boost COMPONENTS python serialization )
find_package( PythonLibs REQUIRED )
if( NOT(Boost_FOUND AND PYTHONLIBS_FOUND AND CPLEX_FOUND) )
  message(FATAL_ERROR "Python wrapper requires boost-python, boost-serialization, the python dev libs, and cplex. Turn off python wrapper if you can't satisfy the dependencies.")
endif()

if(Boost_FOUND AND PYTHONLIBS_FOUND)
  include_directories(${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})
  link_directories( ${Boost_LIBRARY_DIRS} ${PYTHON_LIBRARIES})
endif()

# set definition for vigra random forest
ADD_DEFINITIONS(-DHasHDF5)

##
# wrapper
##

set(PYPGMLINK_SRCS
  pyfield_of_view.cpp
  pyhypotheses.cpp
  pypgmlink.cpp
  pytrack.cpp
  pytraxels.cpp
  )

add_library(pypgmlink SHARED ${PYPGMLINK_SRCS})
target_link_libraries(pypgmlink pgmlink ${CPLEX_LIBRARIES} ${Boost_LIBRARIES} )
set_target_properties(pypgmlink PROPERTIES PREFIX "" OUTPUT_NAME "pgmlink")

##
# install
##

# this will set the install dir to $PYTHON_SITE_PACKAGES
# make use of ${PYTHON_EXECUTABLE} instead of /usr/bin/python

FIND_PACKAGE(PythonInterp)
IF(NOT DEFINED PYPGMLINK_INSTALL_DIR OR PYPGMLINK_INSTALL_DIR MATCHES "^$")
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c 
                    "from distutils.sysconfig import *; print get_python_lib(1)"
                    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
    FILE(TO_CMAKE_PATH ${PYTHON_SITE_PACKAGES} PYPGMLINK_INSTALL_DIR)
ENDIF()
SET(PYPGMLINK_INSTALL_DIR ${PYPGMLINK_INSTALL_DIR}
    CACHE PATH "where to install the pgmlink Python package" FORCE)
# this is the install path relative to CMAKE_INSTALL_PREFIX,
# use this in INSTALL() commands to get packaging right
FILE(RELATIVE_PATH PYPGMLINK_INSTALL_DIR ${CMAKE_INSTALL_PREFIX} ${PYPGMLINK_INSTALL_DIR})

INSTALL(TARGETS pypgmlink
        DESTINATION ${PYPGMLINK_INSTALL_DIR})
