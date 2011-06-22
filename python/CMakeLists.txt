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
