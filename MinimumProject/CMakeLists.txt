cmake_minimum_required(VERSION 2.6)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -std=c++11 -Wall")
set(CMAKE_BUILD_TYPE debug)
project(minimumProject)

add_subdirectory(Library)
add_subdirectory(Tests)