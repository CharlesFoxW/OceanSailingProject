cmake_minimum_required(VERSION 3.7)
project(OceanSailing)

#########################################################
# SET COMPILER
#########################################################
set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -std=c++11)

#########################################################
# FIND GLUT
#########################################################
find_package(GLUT REQUIRED)
include_directories(${GLUT_INCLUDE_DIRS})
link_directories(${GLUT_LIBRARY_DIRS})
add_definitions(${GLUT_DEFINITIONS})
if(NOT GLUT_FOUND)
    message(ERROR " GLUT not found!")
endif(NOT GLUT_FOUND)

#########################################################
# FIND OPENGL
#########################################################
find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})
link_directories(${OpenGL_LIBRARY_DIRS})
add_definitions(${OpenGL_DEFINITIONS})
if(NOT OPENGL_FOUND)
    message(ERROR " OPENGL not found!")
endif(NOT OPENGL_FOUND)

set(SOURCE_FILES main.cpp openGLHeaders.h geometry.h geometry.cpp oceanSailing.h oceanSailing.cpp)

add_executable(OceanSailing ${SOURCE_FILES})
#########################################################
# LINK LIBRARIES
#########################################################
target_link_libraries(OceanSailing ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})