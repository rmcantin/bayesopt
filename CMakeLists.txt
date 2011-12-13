PROJECT(BayesianOptimization CXX)
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

SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -fPIC -fwrapv -fno-strict-aliasing")
SET(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -fPIC -fwrapv -fno-strict-aliasing")

SET(CMAKE_CXX_FLAGS_DEBUG "-g -pg")
SET(CMAKE_C_FLAGS_DEBUG "-g -pg")
  
# enable warnings
ADD_DEFINITIONS("-Wall")

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

OPTION(BUILD_SHARED_LIBS "Build shared libraries?" ON)
option(BUILD_PYTHON_INTERFACE "Build Python Interface?" ON)


if(BUILD_PYTHON_INTERFACE)
  include(FindPythonInterp)
  ## check python
  find_package(PythonLibs) # using PYTHON_INCLUDE_PATH instead of PYTHON_INCLUDE_DIR
  if( NOT PYTHON_EXECUTABLE )
    # look specifically for 2.7
    find_program(PYTHON_EXECUTABLE NAMES python2.7 python27 python PATHS [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.7\\InstallPath])
  endif()

  if( PYTHON_EXECUTABLE )
    # architecture independent
    execute_process(
      COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print get_python_lib(0)"
      OUTPUT_VARIABLE _python_sitepackage OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE _python_failed0)
    # architexture dependent
    execute_process(
      COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print get_python_lib(1)"
      OUTPUT_VARIABLE _python_distpackage OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE _python_failed1)
    execute_process(
      COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; from os.path import relpath; print relpath(get_python_lib(1,prefix='${CMAKE_INSTALL_PREFIX}'),'${CMAKE_INSTALL_PREFIX}')"
      OUTPUT_VARIABLE OPENRAVE_PYTHON_INSTALL_DIR OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE _python_failed2)

    if( ${_python_failed0} EQUAL 0 AND ${_python_failed1} EQUAL 0 AND ${_python_failed2} EQUAL 0 )
      if( NOT IS_DIRECTORY "${_python_distpackage}/numpy/core/include" )
        if( NOT IS_DIRECTORY "${_python_sitepackage}/numpy/core/include" )
          message(STATUS "cannot find python numpy dir!")
          set(PYTHON_EXECUTABLE)
        else()
          set(PYTHON_INCLUDE_PATH ${PYTHON_INCLUDE_PATH} ${_python_sitepackage_arch}/numpy/core/include)
        endif()
      else()
        set(PYTHON_INCLUDE_PATH ${PYTHON_INCLUDE_PATH} ${_python_distpackage}/numpy/core/include)
      endif()

      # get the major.minor python version
      execute_process(
        COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print('%d.%d'%sys.version_info[0:2])"
        OUTPUT_VARIABLE _python_version OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE _python_failed)
      if( ${_python_failed} EQUAL 0 )
        string(REGEX REPLACE "[\r\n]" "" PYTHON_MAJORMINOR_VERSION "${_python_version}")
      else()
        message(STATUS "failed to get python version")
      endif()
    else()
      message(STATUS "failed to get python site-package directories via get_pytho_lib")
      set(PYTHON_EXECUTABLE)
    endif()
  endif()
endif()

IF(BUILD_PYTHON_INTERFACE)
  FIND_PACKAGE(Numpy REQUIRED)
  SET(PYTHON_LIB  ${PYTHON_LIBRARIES} )
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH}
                                       ${PYTHON_NUMPY_INCLUDE_DIR} 
                                       ${PYTHON_INCLUDE_DIR})
ENDIF(BUILD_PYTHON_INTERFACE)

FIND_LIBRARY(NLOPT nlopt)

IF(NLOPT MATCHES NLOPT-NOTFOUND)
  SET(USE_DIRECT_FORTRAN ON CACHE BOOL "Use Fortran Code of DIRECT?")
ELSE(NLOPT MATCHES NLOPT-NOTFOUND)
  SET(USE_DIRECT_FORTRAN OFF CACHE BOOL "Use Fortran Code of DIRECT?")
ENDIF(NLOPT MATCHES NLOPT-NOTFOUND)

SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE(cmake/UseDoxygen)


CONFIGURE_FILE(
  ./include/krig_config.h.in
  ./include/krig_config.h
  )

# with SET() command you can change variables or define new ones
# here we define SAMPLE_SRCS variable that contains a list of all .cpp files
# note that we don't need \ at the end of line
SET( BAYESOPT_SRCS
  ./src/bayesopt.cpp
  ./src/inneroptimization.cpp
  ./src/gaussprocess.cpp
  ./src/basicgaussprocess.cpp
  ./src/studenttprocess.cpp
  ./wrappers/bayesoptwpr.cpp
  )

IF(BUILD_PYTHON_INTERFACE)
  SET(PYTHON_SRC  ./python/bayesopt_py.cpp)
ELSE(BUILD_PYTHON_INTERFACE)
  SET(PYTHON_SRC )
ENDIF(BUILD_PYTHON_INTERFACE)


INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/include 
                     ${CMAKE_SOURCE_DIR}/wrappers
		     ${CMAKE_SOURCE_DIR}/utils )

find_package( Boost REQUIRED )
if(Boost_FOUND)
   include_directories(${Boost_INCLUDE_DIRS})
else(Boost_FOUND)
   find_library(Boost boost PATHS /opt/local/lib)
   include_directories(${Boost_LIBRARY_PATH})
endif()

LINK_DIRECTORIES( ${CMAKE_SOURCE_DIR}/lib )

IF(USE_DIRECT_FORTRAN)
  ADD_SUBDIRECTORY(direct)
  SET( OPT_SOURCE ./wrappers/directwpr.cpp)
ELSE(USE_DIRECT_FORTRAN)
  SET( OPT_SOURCE ./wrappers/nloptwpr.cpp)
ENDIF(USE_DIRECT_FORTRAN)

ADD_LIBRARY(bayesopt
  ${BAYESOPT_SRCS}   
  ${PYTHON_SRC}
  ${OPT_SOURCE}
  )

SET_TARGET_PROPERTIES(bayesopt PROPERTIES PREFIX "" SUFFIX ".so")

IF(USE_DIRECT_FORTRAN)
  SET(EXT_LIBS DIRect)
  add_dependencies(bayesopt DIRect)
ELSE(USE_DIRECT_FORTRAN)
  SET(EXT_LIBS ${NLOPT})
ENDIF(USE_DIRECT_FORTRAN)

TARGET_LINK_LIBRARIES(bayesopt
  ${EXT_LIBS}
  ${PYTHON_LIB}
  )

#Naive test
ADD_EXECUTABLE(bayesopttest ./app/testbayesopt.cpp)
add_dependencies(bayesopttest bayesopt)
TARGET_LINK_LIBRARIES(bayesopttest bayesopt)

#1D test
ADD_EXECUTABLE(bayesoptonedtest ./app/testbopt_oned.cpp)
add_dependencies(bayesoptonedtest bayesopt)
TARGET_LINK_LIBRARIES(bayesoptonedtest bayesopt)

#Branin
ADD_EXECUTABLE(branintest ./app/testbranin.cpp )
add_dependencies(branintest bayesopt)
TARGET_LINK_LIBRARIES(branintest bayesopt )

#Test for random number generator
ADD_EXECUTABLE(randtest ./app/testrand.cpp)

INSTALL(FILES ./include/bayesopt.hpp ./wrappers/bayesoptwpr.h 
  DESTINATION include)

INSTALL(TARGETS bayesopt
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
