add_definitions(-DDATA_PATH="${PROJECT_SOURCE_DIR}/matlabPlots/CPP/")

include_directories("${PROJECT_SOURCE_DIR}/OneDimElliptic")
include_directories("${PROJECT_SOURCE_DIR}/MCMCLib")

add_executable(InversePDE PDE_InverseProblem.cpp)
target_link_libraries(InversePDE MCMC OneDimEllipticSolver)