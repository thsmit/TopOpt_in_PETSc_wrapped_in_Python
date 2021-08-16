# Author: Thijs Smit, Dec 2020
# Copyright (C) 2020 ETH Zurich

# Disclaimer:
# The authors reserves all rights but does not guaranty that the code is
# free from errors. Furthermore, we shall not be liable in any event
# caused by the use of the program.

# check minimum required cmake version
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# set the compiler
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

# set the project name
project(topoptlib VERSION 0.1.0
                  DESCRIPTION "A Python wrapper and extension of the functionality of the TopOpt_in_PETSc framework")

list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake-modules)

find_package(PythonInterp 3.7 REQUIRED)
find_package(PythonLibs 3.7 REQUIRED)
find_package(PETSc REQUIRED)

message (STATUS "PETSC COMPILER ${PETSC_COMPILER}")

add_definitions (-std=c++0x)
include_directories("${CMAKE_SOURCE_DIR}/topoptlib/src/include" "${CMAKE_BINARY_DIR}/include" "${PETSC_INCLUDES}" "${PYTHON_INCLUDE_DIRS}")

add_definitions (${PETSC_DEFINITIONS})
add_definitions (-g)

message (STATUS "PETSC_DEFINITIONS ${PETSC_DEFINITIONS}")

file(GLOB SOURCES "topoptlib/src/*.cc")
add_library(topoptlib SHARED ${SOURCES})

target_link_libraries(topoptlib ${PETSC_LIBRARIES})
target_link_libraries(topoptlib ${PYTHON_LIBRARIES})

set_target_properties(
    topoptlib
    PROPERTIES
        PREFIX ""
        OUTPUT_NAME "topoptlib"
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}
    )
