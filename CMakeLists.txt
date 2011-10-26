PROJECT(BayesianOptimization CXX)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#SET(CMAKE_BUILD_TYPE Release)
#SET(CMAKE_BUILD_TYPE debug)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING 
      "Choose the type of build, options are: Debug Release 
      RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -fPIC -fwrapv -fno-strict-aliasing")
SET(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -fPIC -fwrapv -fno-strict-aliasing")

SET(CMAKE_CXX_FLAGS_DEBUG "-g")
SET(CMAKE_C_FLAGS_DEBUG "-g")
  
# enable warnings
ADD_DEFINITIONS("-Wall -ansi -pedantic")

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
    # look specifically for 2.6
    find_program(PYTHON_EXECUTABLE NAMES python2.6 python PATHS [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.6\\InstallPath])
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

IF(BUILD_PYTHON_INTERFACE)
  SET(PYHTON_LIB  ${PYTHON_LIBRARIES} )
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
ELSE(BUILD_PYTHON_INTERFACE)
  SET(PYHTON_LIB )
ENDIF(BUILD_PYTHON_INTERFACE)

endif()

#FIND_PACKAGE(PythonLibs 2.6)
#IF(PYTHONLIBS_FOUND)
#  SET(BUILD_PYTHON_INTERFACE ON CACHE BOOL "Build Python Interface?")
#  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
#ELSE(PYTHONLIBS_FOUND)
#  SET(BUILD_PYTHON_INTERFACE OFF CACHE BOOL "Build Python Interface?")
#ENDIF(PYTHONLIBS_FOUND)


FIND_LIBRARY(NLOPT nlopt)

IF(NLOPT MATCHES NLOPT-NOTFOUND)
  SET(USE_DIRECT_FORTRAN ON CACHE BOOL "Use Fortran Code of DIRECT?")
ELSE(NLOPT MATCHES NLOPT-NOTFOUND)
  SET(USE_DIRECT_FORTRAN OFF CACHE BOOL "Use Fortran Code of DIRECT?")
ENDIF(NLOPT MATCHES NLOPT-NOTFOUND)

SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE(UseDoxygen)


CONFIGURE_FILE(
  ./include/krig_config.h.in
  ./include/krig_config.h
  )

# with SET() command you can change variables or define new ones
# here we define SAMPLE_SRCS variable that contains a list of all .cpp files
# note that we don't need \ at the end of line
SET( KRIGGING_SRCS
  ./src/krigging.cpp
  ./src/inneroptimization.cpp
  ./src/gaussprocess.cpp
  ./src/basicgaussprocess.cpp
  ./src/studenttprocess.cpp
  ./wrappers/krigwpr.cpp
  )

IF(BUILD_PYTHON_INTERFACE)
  SET(PYHTON_SRC  ./python/bayesopt.cpp)
ELSE(BUILD_PYTHON_INTERFACE)
  SET(PYHTON_SRC )
ENDIF(BUILD_PYTHON_INTERFACE)


INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/include 
                     ${CMAKE_SOURCE_DIR}/wrappers
		     ${CMAKE_SOURCE_DIR}/utils )

find_package( Boost REQUIRED )
   if(Boost_FOUND)
      include_directories(${Boost_INCLUDE_DIRS})
   endif()

LINK_DIRECTORIES( ${CMAKE_SOURCE_DIR}/lib )

IF(USE_DIRECT_FORTRAN)
  ADD_SUBDIRECTORY(direct)
  SET( OPT_SOURCE ./wrappers/directwpr.cpp)
ELSE(USE_DIRECT_FORTRAN)
  SET( OPT_SOURCE ./wrappers/nloptwpr.cpp)
ENDIF(USE_DIRECT_FORTRAN)

ADD_LIBRARY(bayesopt
  ${KRIGGING_SRCS}   
  ${PYHTON_SRC}
  ${OPT_SOURCE}
  )

SET_TARGET_PROPERTIES(bayesopt PROPERTIES PREFIX "")

IF(USE_DIRECT_FORTRAN)
  SET(EXT_LIBS DIRect)
  add_dependencies(bayesopt DIRect)
ELSE(USE_DIRECT_FORTRAN)
  SET(EXT_LIBS ${NLOPT})
ENDIF(USE_DIRECT_FORTRAN)

ADD_EXECUTABLE(krigtest
  ./app/testkrig.cpp
  )

add_dependencies(krigtest bayesopt)

TARGET_LINK_LIBRARIES(bayesopt
  ${EXT_LIBS}
  )

TARGET_LINK_LIBRARIES(krigtest
  bayesopt
  ${PYHTON_LIB}
  )

ADD_EXECUTABLE(randtest ./app/testrand.cpp)

INSTALL(FILES ./include/krigging.hpp ./wrappers/krigwpr.h 
  DESTINATION include)

INSTALL(TARGETS bayesopt
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)