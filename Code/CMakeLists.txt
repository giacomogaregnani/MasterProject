cmake_minimum_required(VERSION 2.6)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopenmp -std=c++11")
set(CMAKE_BUILD_TYPE Debug)
project(projecttest)

include_directories(${CMAKE_SOURCE_DIR}/eigen)

add_subdirectory(Solver)
add_subdirectory(Test)
