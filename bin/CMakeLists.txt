add_executable(bot-trainer bot-trainer.cpp)
target_link_libraries(bot-trainer hdf5 hdf5_hl vigraimpex ilocplex cplex concert bot)