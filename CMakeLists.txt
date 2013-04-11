PROJECT(BayesOpt CXX)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} 
		      ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /opt/local/lib /opt/local/Library)
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} /opt/local/include
/opt/local/Library)
set(CMAKE_PROGRAM_PATH ${CMAKE_PROGRAM_PATH} /opt/local/bin/ /opt/local/Library)


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING 
      "Choose the type of build, options are: Debug Release 
      RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

#SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wdouble-promotion -Wpedantic -Wfloat-equal -Wno-endif-labels -Wshadow")
#SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} -Wall -Wextra -Wdouble-promotion -Wpedantic -Wfloat-equal -Wno-endif-labels -Wshadow")

# SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG ${CMAKE_CXX_FLAGS}")
# SET(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG ${CMAKE_C_FLAGS}")

# SET(CMAKE_CXX_FLAGS_DEBUG "-g -pg ${CMAKE_CXX_FLAGS}")
# SET(CMAKE_C_FLAGS_DEBUG "-g -pg ${CMAKE_C_FLAGS}")
  
# enable warnings
# ADD_DEFINITIONS("-Wall")
# ADD_DEFINITIONS("-Wextra -Wdouble-promotion -pedantic -Wfloat-equal -Wno-endif-labels -Wshadow")

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

option(BUILD_PYTHON_INTERFACE "Build Python Interface?" OFF)
option(MATLAB_COMPATIBLE "Build library compatible with Matlab?" ON)

if(BUILD_PYTHON_INTERFACE)
  INCLUDE(PythonMagic)
  SET(PYTHON_LIB  ${PYTHON_LIBRARIES} )
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
  SET(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries?" FORCE)
ELSE()
  SET(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries?")
ENDIF(BUILD_PYTHON_INTERFACE)

find_package( Boost REQUIRED )
if(Boost_FOUND)
   include_directories(${Boost_INCLUDE_DIRS})
else(Boost_FOUND)
   find_library(Boost boost PATHS /opt/local/lib)
   include_directories(${Boost_LIBRARY_PATH})
endif()

FIND_LIBRARY(NLOPT nlopt)

IF(NLOPT MATCHES NLOPT-NOTFOUND)
  SET(BUILD_NLOPT ON CACHE BOOL "Build included version of NLOPT?")
ELSE(NLOPT MATCHES NLOPT-NOTFOUND)
  SET(BUILD_NLOPT OFF CACHE BOOL "Build included version of NLOPT?")
ENDIF(NLOPT MATCHES NLOPT-NOTFOUND)

#SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE(UseDoxygen)

SET( BAYESOPT_SRCS
  ./src/bayesoptcont.cpp
  ./src/bayesoptdisc.cpp
  ./src/bayesoptbase.cpp
  ./src/inneroptimization.cpp
  ./src/nonparametricprocess.cpp
  ./src/hierarchical_gaussian_process.cpp
  ./src/gaussian_process.cpp
  ./src/gaussian_process_ml.cpp
#  ./src/gaussian_process_ign.cpp
#  ./src/student_t_process.cpp
  ./src/parameters.cpp
  ./src/kernel_functors.cpp
  ./src/criteria_functors.cpp
  ./src/criteria_combined.cpp
  ./src/mean_functors.cpp
  ./src/gauss_distribution.cpp
  ./src/student_t_distribution.cpp
  )

SET(UTILS_SRC
  ./utils/parser.cpp
  )

SET(WRAPPPERS_SRC 
  ./wrappers/nloptwpr.cpp
  ./wrappers/bayesoptwpr.cpp
  )

IF(BUILD_PYTHON_INTERFACE)
  SET(PYTHON_SRC  ./python/bayesopt.cpp)
ELSE(BUILD_PYTHON_INTERFACE)
  SET(PYTHON_SRC )
ENDIF(BUILD_PYTHON_INTERFACE)


INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/include 
                     ${CMAKE_SOURCE_DIR}/wrappers
		     ${CMAKE_SOURCE_DIR}/utils )


LINK_DIRECTORIES( ${CMAKE_SOURCE_DIR}/lib )

# IF(BUILD_NLOPT)
#   ADD_SUBDIRECTORY(nlopt)
#   include_directories(${CMAKE_SOURCE_DIR}/nlopt/api)
# ENDIF(BUILD_NLOPT)

IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(bayesopt SHARED ${BAYESOPT_SRCS}   
    ${WRAPPPERS_SRC} ${UTILS_SRC} ${PYTHON_SRC} )
  IF(WIN32)
    ADD_DEFINITIONS(-DBAYESOPT_DLL)	
    # In new versions of CMAKE they use a different system and the
    # symbol is not defined
    ADD_DEFINITIONS(-Dbayesopt_EXPORT )
  ELSE()
    SET_TARGET_PROPERTIES(bayesopt PROPERTIES PREFIX "" SUFFIX ".so")
  ENDIF()
ELSE()
  ADD_LIBRARY(bayesopt STATIC ${BAYESOPT_SRCS}   
    ${WRAPPPERS_SRC} ${UTILS_SRC} ${PYTHON_SRC} )
ENDIF()

IF((BUILD_SHARED_LIBS OR MATLAB_COMPATIBLE) AND NOT WIN32)
  SET_TARGET_PROPERTIES(bayesopt PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF()
  

IF(BUILD_NLOPT)
  ADD_SUBDIRECTORY(nlopt)
  include_directories(${CMAKE_SOURCE_DIR}/nlopt/api)
  add_dependencies(bayesopt nlopt)
  SET(EXT_LIBS nlopt)
ELSE(BUILD_NLOPT)
  SET(EXT_LIBS ${NLOPT})
ENDIF(BUILD_NLOPT)

TARGET_LINK_LIBRARIES(bayesopt
  ${EXT_LIBS} ${PYTHON_LIB} )

#Continuous test
ADD_EXECUTABLE(bo_cont ./app/bo_cont.cpp)
add_dependencies(bo_cont bayesopt)
TARGET_LINK_LIBRARIES(bo_cont bayesopt)

#Discrete test
ADD_EXECUTABLE(bo_disc ./app/bo_disc.cpp)
add_dependencies(bo_disc bayesopt)
TARGET_LINK_LIBRARIES(bo_disc bayesopt)

#1D test
ADD_EXECUTABLE(bo_oned ./app/bo_oned.cpp)
add_dependencies(bo_oned bayesopt)
TARGET_LINK_LIBRARIES(bo_oned bayesopt)

#Branin
ADD_EXECUTABLE(bo_branin ./app/bo_branin.cpp )
add_dependencies(bo_branin bayesopt)
TARGET_LINK_LIBRARIES(bo_branin bayesopt)

#Test for random number generator
#ADD_EXECUTABLE(randtest ./app/testrand.cpp)

#Test for parsers
#ADD_EXECUTABLE(parsetest ./utils/parser.cpp ./app/testparser.cpp)


INSTALL(FILES 
  ./include/bayesoptcont.hpp 
  ./include/bayesoptdisc.hpp 
  ./wrappers/bayesoptwpr.h 
  DESTINATION include
)

INSTALL(
  TARGETS bayesopt
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
