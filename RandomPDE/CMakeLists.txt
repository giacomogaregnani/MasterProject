cmake_minimum_required(VERSION 2.6)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -std=c++11 -Wall")
set(CMAKE_BUILD_TYPE Release)
project(randomPDE)

include_directories(${CMAKE_SOURCE_DIR}/eigen)

add_subdirectory(FirstTests)
