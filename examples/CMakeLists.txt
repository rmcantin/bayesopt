CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#Continuous test
ADD_EXECUTABLE(bo_cont ./bo_cont.cpp)
add_dependencies(bo_cont bayesopt)
TARGET_LINK_LIBRARIES(bo_cont bayesopt)

#Discrete test
ADD_EXECUTABLE(bo_disc ./bo_disc.cpp)
add_dependencies(bo_disc bayesopt)
TARGET_LINK_LIBRARIES(bo_disc bayesopt)

#1D test
ADD_EXECUTABLE(bo_oned ./bo_oned.cpp)
add_dependencies(bo_oned bayesopt)
TARGET_LINK_LIBRARIES(bo_oned bayesopt)

#Branin
ADD_EXECUTABLE(bo_branin ./bo_branin.cpp )
add_dependencies(bo_branin bayesopt)
TARGET_LINK_LIBRARIES(bo_branin bayesopt)

enable_language(C)  #For OpenGL and bo_display demo

#Display test
find_package(GLUT)
find_package(OpenGL)
if(OPENGL_FOUND AND GLUT_FOUND)
  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/matplotpp)

  include_directories(${GLUT_INCLUDE_DIRS})
  link_directories(${GLUT_LIBRARY_DIRS})
  add_definitions(${GLUT_DEFINITIONS})

  include_directories(${OpenGL_INCLUDE_DIRS})
  link_directories(${OpenGL_LIBRARY_DIRS})
  add_definitions(${OpenGL_DEFINITIONS})
 
  ADD_EXECUTABLE(bo_display 
    ./bo_display.cpp 
    ../matplotpp/matplotpp.cc 
    ../matplotpp/gl2ps.c)
  add_dependencies(bo_display bayesopt)
  TARGET_LINK_LIBRARIES(bo_display bayesopt ${GLUT_LIBRARY} ${OPENGL_LIBRARY})
endif()
