cmake_minimum_required(VERSION 2.6)

PROJECT(SuperPoints)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(SuperPoints SuperPoints.cxx vtkSuperPoints.cxx)
TARGET_LINK_LIBRARIES(SuperPoints vtkHybrid vtkInfovis)
