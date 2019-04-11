cmake_minimum_required(VERSION 2.6)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -std=c++11 -Wall")
set(CMAKE_BUILD_TYPE Release)
project(randomPDE)

include_directories(${CMAKE_SOURCE_DIR}/eigen)
include_directories(${CMAKE_SOURCE_DIR}/getpot-c++)

add_subdirectory(MCMC)
add_subdirectory(NumericalSolver)
add_subdirectory(SDETests)
add_subdirectory(InverseProblemsTests)
add_subdirectory(PFTEsts)
add_subdirectory(ParticleFilter)