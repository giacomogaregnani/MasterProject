add_definitions(-DDATA_PATH="${PROJECT_SOURCE_DIR}/matlabPlots/RandomFEM_BAYESIAN/")

include_directories("${PROJECT_SOURCE_DIR}/RKSolver")
include_directories("${PROJECT_SOURCE_DIR}/RKRandomStep")
include_directories("${PROJECT_SOURCE_DIR}/OneDimElliptic")
include_directories("${PROJECT_SOURCE_DIR}/MCMCLib")

add_executable(InversePDEMulti PDE_InverseProblem_InfDim.cpp)
target_link_libraries(InversePDEMulti MCMC OneDimEllipticSolver)

add_executable(InversePDEFinite PDE_InverseProblem_FinDim.cpp)
target_link_libraries(InversePDEFinite MCMC OneDimEllipticSolver)
