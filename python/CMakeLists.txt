cmake_minimum_required(VERSION 2.8)
message( "\nConfiguring python wrapper:" )

# dependencies
find_package( Boost REQUIRED COMPONENTS python serialization )
if(WIN32)
  set(PYTHON_EXECUTABLE "" CACHE PATH "")
  set(PYTHON_INCLUDE_DIRS "" CACHE PATH "")
  set(PYTHON_LIBRARIES "" CACHE PATH "")
else()
  find_package( PythonLibs REQUIRED )
  find_package( PythonInterp REQUIRED )
endif()

include_directories( ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS}/python2.7 )

# pypgmlink
set(PYPGMLINK_SRCS
  pyfield_of_view.cpp
  pyhypotheses.cpp
  pypgmlink.cpp
  pytrack.cpp
  pytraxels.cpp
  )

add_library(pypgmlink SHARED ${PYPGMLINK_SRCS})#
target_link_libraries(pypgmlink pgmlink ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} )
if(WIN32)
  set_target_properties(pypgmlink PROPERTIES PREFIX "" OUTPUT_NAME "pgmlink" SUFFIX ".pyd")
else()
  set_target_properties(pypgmlink PROPERTIES PREFIX "" OUTPUT_NAME "pgmlink" SUFFIX ".so")
endif()

# install

# this will set the install dir to $PYTHON_SITE_PACKAGES
# make use of ${PYTHON_EXECUTABLE} instead of /usr/bin/python
if(NOT DEFINED PYPGMLINK_INSTALL_DIR OR PYPGMLINK_INSTALL_DIR MATCHES "^$")
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c 
                    "from distutils.sysconfig import *; print get_python_lib(1)"
                    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
    file(TO_CMAKE_PATH ${PYTHON_SITE_PACKAGES} PYPGMLINK_INSTALL_DIR)
endif()
set(PYPGMLINK_INSTALL_DIR ${PYPGMLINK_INSTALL_DIR}
    CACHE PATH "where to install the pgmlink Python package" FORCE)

# this is the install path relative to CMAKE_INSTALL_PREFIX,
# use this in INSTALL() commands to get packaging right
file(RELATIVE_PATH PYPGMLINK_INSTALL_DIR ${CMAKE_INSTALL_PREFIX} ${PYPGMLINK_INSTALL_DIR})

install(TARGETS pypgmlink
	RUNTIME DESTINATION ${PYPGMLINK_INSTALL_DIR}
        LIBRARY DESTINATION ${PYPGMLINK_INSTALL_DIR})
