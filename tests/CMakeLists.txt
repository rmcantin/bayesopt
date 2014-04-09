CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#Test for random number generator
ADD_EXECUTABLE(randtest ./testrand.cpp)

#Test for parsers
ADD_EXECUTABLE(parsetest ../utils/parser.cpp ./testparser.cpp)

#Test for grid
#ADD_EXECUTABLE(gridtest ./testgrid.cpp)

#Test for MCMC
# ADD_EXECUTABLE(mcmctest ../src/mcmc.cpp ./testmcmc.cpp)

#Continuous test
# ADD_EXECUTABLE(bo_cont ./bo_cont.cpp)
# add_dependencies(bo_cont bayesopt)
# TARGET_LINK_LIBRARIES(bo_cont bayesopt)
