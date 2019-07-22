cmake_minimum_required(VERSION 2.8)

SET(CMAKE_BUILD_TYPE "Debug")
#SET(VTK_DIR "/home/tah09e/code/workspace/vtkNightly/build")
#SET(VTK_DIR "/global/homes/t/tah09e/code/workspace/vtk5.10.1/build")
SET(VTK_DIR "/home/tah09e/code/workspace/vtk5.10.0/build")

PROJECT(LAVAflow) 

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc) 
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/Tests) 
