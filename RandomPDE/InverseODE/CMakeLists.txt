add_definitions(-DDATA_PATH="${PROJECT_SOURCE_DIR}/matlabPlots/MCMC_ODE/")

include_directories("${PROJECT_SOURCE_DIR}/RKSolver")
include_directories("${PROJECT_SOURCE_DIR}/RKRandomStep")
include_directories("${PROJECT_SOURCE_DIR}/OneDimElliptic")
include_directories("${PROJECT_SOURCE_DIR}/MCMCLib")
include_directories("${PROJECT_SOURCE_DIR}/SMCLib")

add_executable(InverseODE ODE_InverseProblem.cpp ReadObservations.hpp ReadObservations.cpp)
target_link_libraries(InverseODE MCMC SMC RungeKuttaSolver RandomStep)

add_executable(InverseODEFixed ODE_InverseProblemFixed.cpp ReadObservations.hpp ReadObservations.cpp)
target_link_libraries(InverseODEFixed MCMC SMC RungeKuttaSolver RandomStep)

add_executable(InverseODEIC ODE_InverseProblemIC.cpp ReadObservations.hpp ReadObservations.cpp)
target_link_libraries(InverseODEIC MCMC RungeKuttaSolver RandomStep)

add_executable(InverseODEICEE ODE_InverseProblemICErrEst.cpp ReadObservations.hpp ReadObservations.cpp)
target_link_libraries(InverseODEICEE MCMC RungeKuttaSolver RandomStep)

#add_executable(InfInverseODE ODE_InfInverseProblem.cpp ReadObservations.hpp ReadObservations.cpp)
#target_link_libraries(InfInverseODE MCMC SMC RungeKuttaSolver RandomStep)

add_executable(GenObs GenerateObservations.cpp)
target_link_libraries(GenObs RungeKuttaSolver)

add_executable(ODE_SMC ODE_SMC.cpp ReadObservations.hpp ReadObservations.cpp)
target_link_libraries(ODE_SMC SMC RungeKuttaSolver RandomStep)

add_executable(KL TestKL.cpp)
target_link_libraries(KL MCMC)
