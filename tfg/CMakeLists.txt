CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#Optimization with system processes
ADD_EXECUTABLE(system_calls_branin ./system_calls_branin.cpp)
add_dependencies(system_calls_branin bayesopt)
TARGET_LINK_LIBRARIES(system_calls_branin bayesopt)
