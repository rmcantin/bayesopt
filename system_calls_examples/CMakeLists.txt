CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#Optimization with system calls
ADD_EXECUTABLE(branin_system_calls ./branin_system_calls.cpp)
add_dependencies(branin_system_calls bayesopt)
TARGET_LINK_LIBRARIES(branin_system_calls bayesopt)

#Optimization with xml and system calls
ADD_EXECUTABLE(branin_xml ./branin_xml.cpp)
add_dependencies(branin_xml bayesopt)
TARGET_LINK_LIBRARIES(branin_xml bayesopt)
