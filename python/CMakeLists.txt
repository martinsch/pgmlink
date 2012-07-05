#
# (c) Bernhard X. Kausler, 2010
#

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake_extensions/)

####################
### dependencies ###
####################

find_package( VIGRA REQUIRED )
find_package( Boost REQUIRED COMPONENTS python )
find_package( PythonLibs REQUIRED )

if(VIGRA_FOUND AND Boost_FOUND AND PYTHONLIBS_FOUND)
  include_directories(${VIGRA_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})
  link_directories( ${Boost_LIBRARY_DIRS} ${VIGRA_IMPEX_LIBRARY_DIR} ${PYTHON_LIBRARIES})
endif()

# assume by default that HDF5 is installed. Set definition for vigra random forest
ADD_DEFINITIONS(-DHasHDF5)



###############
### wrapper ###
###############

set(TRACKING_SRCS ctracking.cpp)
link_directories(${PROJECT_SOURCE_DIR}/../3rd_party/lib/)

add_library(ctracking SHARED ${TRACKING_SRCS})
target_link_libraries(ctracking ${Boost_PYTHON_LIBRARIES} track_lib ilocplex vigraimpex)
set_target_properties(ctracking PROPERTIES PREFIX "" OUTPUT_NAME "ctracking")


# install




# this will set ctracking install dir to $PYTHON_SITE_PACKAGES
# make use of ${PYTHON_EXECUTABLE} instead of /usr/bin/python


FIND_PACKAGE(PythonInterp)
IF(NOT DEFINED CTRACKING_INSTALL_DIR OR CTRACKING_INSTALL_DIR MATCHES "^$")
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c 
                    "from distutils.sysconfig import *; print get_python_lib(1)"
                    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
    FILE(TO_CMAKE_PATH ${PYTHON_SITE_PACKAGES} CTRACKING_INSTALL_DIR)
ENDIF()
SET(CTRACKING_INSTALL_DIR ${CTRACKING_INSTALL_DIR}
    CACHE PATH "where to install the VIGRA Python package" FORCE)
# this is the install path relative to CMAKE_INSTALL_PREFIX,
# use this in INSTALL() commands to get packaging right
FILE(RELATIVE_PATH CTRACKING_INSTALL_DIR ${CMAKE_INSTALL_PREFIX} ${CTRACKING_INSTALL_DIR})


INSTALL(TARGETS ctracking
        DESTINATION ${CTRACKING_INSTALL_DIR}/ctracking)