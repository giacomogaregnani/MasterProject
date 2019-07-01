add_definitions(-DDATA_PATH="${PROJECT_SOURCE_DIR}/matlabPlots/testParticle/")

include_directories("${PROJECT_SOURCE_DIR}/NumericalSolver")
include_directories("${PROJECT_SOURCE_DIR}/ParticleFilter")

add_executable(testPF
               generateObservations.hpp generateObservations.cpp
               PFTest.cpp)
target_link_libraries(testPF ParFil SDESolver)