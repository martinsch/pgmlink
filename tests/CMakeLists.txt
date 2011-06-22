cmake_minimum_required(VERSION 2.8)

find_package( Boost REQUIRED COMPONENTS unit_test_framework )

if(Boost_FOUND)
  add_definitions(-DBOOST_TEST_DYN_LINK)

  include_directories(
    ${Boost_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/include/
    ${PROJECT_SOURCE_DIR}/../../3rd_party/include
  )

  LINK_DIRECTORIES( 
    ${PROJECT_SOURCE_DIR}/../../3rd_party/lib
  )


  add_executable(ilp_construction_test ilp_construction_test.cpp)
  target_link_libraries(ilp_construction_test track_lib ${Boost_LIBRARIES} ann)

  add_executable(track_test track_test.cpp)
  target_link_libraries(track_test track_lib ${Boost_LIBRARIES} vigraimpex ilocplex)

  add_executable(traxels_test traxels_test.cpp)
  target_link_libraries(traxels_test track_lib ${Boost_LIBRARIES})

  add_executable(energy_test energy_test.cpp)
  target_link_libraries(energy_test track_lib ${Boost_LIBRARIES})

  add_executable(graph_test graph_test.cpp)
  target_link_libraries(graph_test track_lib ${Boost_LIBRARIES})

  add_executable(graphical_model_test graphical_model_test.cpp) 
  target_link_libraries(graphical_model_test track_lib ${Boost_LIBRARIES} vigraimpex ilocplex cplex concert)

  add_executable(mrf_reasoner_test mrf_reasoner_test.cpp) 
  target_link_libraries(mrf_reasoner_test track_lib ${Boost_LIBRARIES} ilocplex cplex concert)

  add_executable(hypotheses_test hypotheses_test.cpp) 
  target_link_libraries(hypotheses_test track_lib ${Boost_LIBRARIES} ilocplex cplex concert)

  # assume by default that HDF5 is installed. Set definition for vigra random forest
  ADD_DEFINITIONS(-DHasHDF5)

  add_executable(randomforest_test randomforest_test.cpp)
  target_link_libraries(randomforest_test track_lib ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} hdf5 hdf5_hl vigraimpex)

  add_test(hypotheses hypotheses_test)
  add_test(mrf_reasoner mrf_reasoner_test)
  add_test(graphical_model graphical_model_test)
  add_test(graph graph_test)
  add_test(ilp_construction ilp_construction_test)
  add_test(traxels traxels_test)
  add_test(energy energy_test)
  add_test(track track_test)
  add_test(randomforest randomforest_test)
endif()

