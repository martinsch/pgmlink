cmake_minimum_required(VERSION 2.8)

##
# dependencies
##
find_package( Boost COMPONENTS unit_test_framework )
if( NOT(Boost_FOUND AND HDF5_FOUND AND VIGRA_FOUND AND CPLEX_FOUND) )
  message(STATUS "Boost_FOUND: ${Boost_FOUND}")
  message(STATUS "HDF5_FOUND: ${HDF5_FOUND}")
  message(STATUS "VIGRA_FOUND: ${VIGRA_FOUND}")
  message(STATUS "Cplex_FOUND: ${Cplex_FOUND}")
  message(FATAL_ERROR "boost-test, hdf5, vigra, and cplex required for test but not found. Turn off tests if you can't satisfy the dependencies.")
endif()


add_definitions(-DBOOST_TEST_DYN_LINK)

include_directories(
  ${Boost_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/include/
  )

## autodiscover test sources and add tests
file(GLOB TEST_SRCS *.cpp)
foreach(test_src ${TEST_SRCS})
  if (${test_src} MATCHES "randomforest_test.cpp")
  elseif(${test_src} MATCHES "learning_systemtest.cpp")
  else()
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable( ${test_name} ${test_src} )
    target_link_libraries( ${test_name} pgmlink ${CPLEX_LIBRARIES} ${Boost_LIBRARIES} ${HDF5_LIBRARIES})   
    add_test( ${test_name} ${test_name})
  endif()
endforeach(test_src)

## randomforest lives in the project dir and we add it manually
# assume by default that HDF5 is installed. Set definition for vigra random forest
ADD_DEFINITIONS(-DHasHDF5)
configure_file(randomforest_test.cpp.cmake randomforest_test.cpp)
add_executable(randomforest_test ${PROJECT_BINARY_DIR}/tests/randomforest_test.cpp)

target_link_libraries(randomforest_test pgmlink ${Boost_LIBRARIES} ${CPLEX_LIBRARIES} ${VIGRA_LIBRARIES} ${HDF5_LIBRARIES} )
add_test(randomforest randomforest_test)

## system tests
configure_file( system/learning_systemtest.cpp.cmake learning_systemtest.cpp)
add_executable( learning_systemtest ${PROJECT_BINARY_DIR}/tests/learning_systemtest.cpp )
target_link_libraries( learning_systemtest pgmlink ${CPLEX_LIBRARIES} ${Boost_LIBRARIES} )   
#add_test( learning_systemtest learning_systemtest )

## opengm extensions tests
file(GLOB OGM_TEST_SRCS ext_opengm/*.cxx)
foreach(test_src ${OGM_TEST_SRCS})
  get_filename_component(test_name ${test_src} NAME_WE)
  set(test_name ext_opengm_${test_name})
  add_executable( ${test_name} ${test_src} )
  target_link_libraries( ${test_name} ${CPLEX_LIBRARIES} )
  add_test( ${test_name} ${test_name} )
endforeach(test_src)



