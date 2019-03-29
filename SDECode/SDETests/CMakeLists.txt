add_definitions(-DDATA_PATH="${PROJECT_SOURCE_DIR}/matlabPlots/testSolver/")

include_directories("${PROJECT_SOURCE_DIR}/NumericalSolver")

add_executable(testSolver testSolver.cpp)
target_link_libraries(testSolver SDESolver)

add_executable(testSolverGN testSolverGivenNoise.cpp)
target_link_libraries(testSolverGN SDESolver)
