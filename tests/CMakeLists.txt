CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#Test for random number generator
ADD_EXECUTABLE(randtest ./testrand.cpp)

#Test for parsers
ADD_EXECUTABLE(parsetest ../utils/parser.cpp ./testparser.cpp)

#Test for grid
ADD_EXECUTABLE(gridtest ./testgrid.cpp)

#Test for MCMC
ADD_EXECUTABLE(mcmctest ../src/mcmc_sampler.cpp ./testmcmc.cpp)

#Test for Cholesky decomposition
ADD_EXECUTABLE(choltest ./testchol.cpp)
