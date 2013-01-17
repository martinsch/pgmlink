cmake_minimum_required(VERSION 2.8)

find_package( Boost REQUIRED COMPONENTS unit_test_framework )

if(Boost_FOUND)
  add_definitions(-DBOOST_TEST_DYN_LINK)

  include_directories(
    ${Boost_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/include/
  )

  ## autodiscover test sources and add tests
  file(GLOB TEST_SRCS *.cpp)
  foreach(test_src ${TEST_SRCS})
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable( ${test_name} ${test_src} )
    target_link_libraries( ${test_name} pgmlink ${Boost_LIBRARIES})   
    add_test( ${test_name} ${test_name})
  endforeach(test_src)

  ## randomforest lives in the project dir and we add it manually
  # assume by default that HDF5 is installed. Set definition for vigra random forest
  ADD_DEFINITIONS(-DHasHDF5)
  configure_file(randomforest_test.cpp.cmake randomforest_test.cpp)
  add_executable(randomforest_test ${PROJECT_BINARY_DIR}/tests/randomforest_test.cpp)
  target_link_libraries(randomforest_test pgmlink ${Boost_LIBRARIES} hdf5 hdf5_hl vigraimpex)
  add_test(randomforest randomforest_test)

endif()

