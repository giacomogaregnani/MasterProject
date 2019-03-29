cmake_minimum_required(VERSION 2.6)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -std=c++11 -Wall")
set(CMAKE_BUILD_TYPE Release)
project(randomODE-MasterProject)

include_directories(${CMAKE_SOURCE_DIR}/eigen)
include_directories(${CMAKE_SOURCE_DIR}/getpot-c++)

add_subdirectory(RKSolver)
add_subdirectory(Test)
add_subdirectory(RKRandomStep)
