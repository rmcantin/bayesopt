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

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

option(BAYESOPT_BUILD_EXAMPLES "Build examples and demos?" ON)
option(BAYESOPT_PYTHON_INTERFACE "Build Python interface?" OFF)
option(BAYESOPT_MATLAB_COMPATIBLE "Build library compatible with Matlab?" ON)
option(BAYESOPT_BUILD_SOBOL "Build support for Sobol sequences?" ON)

if(BAYESOPT_PYTHON_INTERFACE)
  INCLUDE(PythonMagic)
  SET(PYTHON_LIB  ${PYTHON_LIBRARIES} )
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
  SET(BAYESOPT_BUILD_SHARED ON CACHE BOOL "Build shared libraries?" FORCE)
ELSE()
  SET(BAYESOPT_BUILD_SHARED OFF CACHE BOOL "Build shared libraries?")
ENDIF(BAYESOPT_PYTHON_INTERFACE)

find_package( Boost REQUIRED )
if(Boost_FOUND)
   include_directories(${Boost_INCLUDE_DIRS})
else(Boost_FOUND)
   find_library(Boost boost PATHS /opt/local/lib)
   include_directories(${Boost_LIBRARY_PATH})
endif()

FIND_LIBRARY(NLOPT nlopt)

IF(NLOPT MATCHES NLOPT-NOTFOUND)
  SET(NLOPT_BUILD ON CACHE BOOL "Build included version of NLOPT?")
ELSE(NLOPT MATCHES NLOPT-NOTFOUND)
  SET(NLOPT_BUILD OFF CACHE BOOL "Build included version of NLOPT?")
ENDIF(NLOPT MATCHES NLOPT-NOTFOUND)

INCLUDE(UseDoxygen)


IF(BAYESOPT_BUILD_SOBOL)
  ADD_DEFINITIONS(-DUSE_SOBOL)	
  SET(SOBOL_SRC
    ./sobol/sobol.cpp
    ./sobol/sobol_i4.cpp
    ./sobol/sobol_i8.cpp
    )
ELSE(BAYESOPT_BUILD_SOBOL)
  SET(SOBOL_SRC )
ENDIF(BAYESOPT_BUILD_SOBOL)


SET( BAYESOPT_SRCS
  ./src/bayesoptcont.cpp
  ./src/bayesoptdisc.cpp
  ./src/bayesoptbase.cpp
  ./src/inneroptimization.cpp
  ./src/nonparametricprocess.cpp
  ./src/hierarchical_gaussian_process.cpp
  ./src/gaussian_process.cpp
  ./src/gaussian_process_ml.cpp
  ./src/gaussian_process_normal.cpp
  ./src/student_t_process_jef.cpp
  ./src/student_t_process_nig.cpp
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
  ./utils/ublas_extra.cpp
  ${SOBOL_SRC}
  )

SET(WRAPPPERS_SRC 
  ./wrappers/nloptwpr.cpp
  ./wrappers/bayesoptwpr.cpp
  )

IF(BAYESOPT_PYTHON_INTERFACE)
  SET(PYTHON_SRC  ./python/bayesopt.cpp)
ELSE(BAYESOPT_PYTHON_INTERFACE)
  SET(PYTHON_SRC)
ENDIF(BAYESOPT_PYTHON_INTERFACE)


INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/include 
                     ${CMAKE_SOURCE_DIR}/wrappers
		     ${CMAKE_SOURCE_DIR}/utils
		     ${CMAKE_SOURCE_DIR}/sobol )


LINK_DIRECTORIES( ${CMAKE_SOURCE_DIR}/lib )


IF(BAYESOPT_BUILD_SHARED)
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

IF((BAYESOPT_BUILD_SHARED OR BAYESOPT_MATLAB_COMPATIBLE) AND NOT WIN32)
  SET_TARGET_PROPERTIES(bayesopt PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF()
  

IF(NLOPT_BUILD)
  ADD_SUBDIRECTORY(nlopt)
  include_directories(${CMAKE_SOURCE_DIR}/nlopt/api)
  add_dependencies(bayesopt nlopt)
  SET(EXT_LIBS nlopt)
ELSE(NLOPT_BUILD)
  SET(EXT_LIBS ${NLOPT})
ENDIF(NLOPT_BUILD)

TARGET_LINK_LIBRARIES(bayesopt
  ${EXT_LIBS} ${PYTHON_LIB} )


IF(BAYESOPT_BUILD_EXAMPLES)
  ADD_SUBDIRECTORY(examples)
endif(BAYESOPT_BUILD_EXAMPLES)

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
