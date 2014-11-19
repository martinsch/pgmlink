cmake_minimum_required(VERSION 2.8)
message( "\nConfiguring tests:" )

# dependencies
find_package( Boost REQUIRED COMPONENTS unit_test_framework )
#find_package( Cplex REQUIRED )

include_directories(
  ${Boost_INCLUDE_DIRS}
  ${OPTIMIZER_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/include/
)

# enable dynamic linking with boost test
add_definitions(-DBOOST_TEST_DYN_LINK)

# autodiscover test sources and add tests
file(GLOB TEST_SRCS *.cpp)
foreach(test_src ${TEST_SRCS})
  # these tests have to be configured first
  if (${test_src} MATCHES "randomforest_test.cpp") 
  elseif(${test_src} MATCHES "learning_systemtest.cpp")
  else()
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable( ${test_name} ${test_src} )
    target_link_libraries( ${test_name} pgmlink ${Boost_LIBRARIES})
    add_test( ${test_name} ${test_name})
  endif()
endforeach(test_src)

# configure random forest test
## assume by default that HDF5 is installed; set definition for vigra random forest
ADD_DEFINITIONS(-DHasHDF5)
configure_file(randomforest_test.cpp.cmake randomforest_test.cpp)
add_executable(randomforest_test ${PROJECT_BINARY_DIR}/tests/randomforest_test.cpp)

target_link_libraries(randomforest_test pgmlink ${Boost_LIBRARIES})
add_test(randomforest randomforest_test)

# system tests
configure_file( system/learning_systemtest.cpp.cmake learning_systemtest.cpp)
add_executable( learning_systemtest ${PROJECT_BINARY_DIR}/tests/learning_systemtest.cpp )
target_link_libraries( learning_systemtest pgmlink ${Boost_LIBRARIES} )
#add_test( learning_systemtest learning_systemtest )

# opengm extensions tests
file(GLOB OGM_TEST_SRCS ext_opengm/*.cxx)
foreach(test_src ${OGM_TEST_SRCS})
  get_filename_component(test_name ${test_src} NAME_WE)
  set(test_name ext_opengm_${test_name})
  add_executable( ${test_name} ${test_src} )
  target_link_libraries( ${test_name} ${OPTIMIZER_LIBRARIES} )
  add_test( ${test_name} ${test_name} )
endforeach(test_src)
