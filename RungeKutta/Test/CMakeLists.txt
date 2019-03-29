add_definitions(-DDATA_PATH="${PROJECT_SOURCE_DIR}/matlabPlots/")

include_directories("${PROJECT_SOURCE_DIR}/RKSolver")
include_directories("${PROJECT_SOURCE_DIR}/RKRandomStep")
include_directories("${PROJECT_SOURCE_DIR}/Utilities")
include_directories("${PROJECT_SOURCE_DIR}/MCMCLib")

add_executable(randomTest RandomTest.cpp)
target_link_libraries(randomTest RandomStep)

add_executable(randomSol RandomSolutions.cpp)
target_link_libraries(randomSol RandomStep)

add_executable(weakConv WeakConvergenceTest.cpp)
target_link_libraries(weakConv RandomStep)

add_executable(strongConv StrongConvergenceTest.cpp)
target_link_libraries(strongConv RandomStep UtilitiesGG)

add_executable(strongConvOne StrongConvergenceOneStep.cpp)
target_link_libraries(strongConvOne RandomStep UtilitiesGG)

add_executable(MCConv MonteCarloTest.cpp)
target_link_libraries(MCConv RandomStep)

add_executable(refSol ReferenceSolution.cpp)
target_link_libraries(refSol RungeKuttaSolver)

add_executable(genObs GenerateObservations.cpp)
target_link_libraries(genObs RungeKuttaSolver)

add_executable(genObsHeat GenerateObservationsHeat.cpp)
target_link_libraries(genObsHeat RungeKuttaSolver)

add_executable(strongConvAdd StrongConvergenceAdd.cpp)
target_link_libraries(strongConvAdd RandomStep UtilitiesGG)

add_executable(weakConvAdd WeakConvergenceAdd.cpp)
target_link_libraries(weakConvAdd RandomStep UtilitiesGG)

add_executable(RandVSAdd RandStepVsAddNoise.cpp)
target_link_libraries(RandVSAdd RandomStep UtilitiesGG)

add_executable(LorenzIC LorenzVariableIC.cpp)
target_link_libraries(LorenzIC RandomStep UtilitiesGG)

add_executable(SpaceMeans SpaceMeans.cpp)
target_link_libraries(SpaceMeans RandomStep UtilitiesGG)

add_executable(TimeMeans TimeMeans.cpp)
target_link_libraries(TimeMeans RandomStep UtilitiesGG)

add_executable(Invariants InvariantCons.cpp)
target_link_libraries(Invariants RandomStep UtilitiesGG)

add_executable(InvariantsExplosion InvariantExplTime.cpp)
target_link_libraries(InvariantsExplosion RandomStep UtilitiesGG)

add_executable(HamError HamiltonianError.cpp)
target_link_libraries(HamError RandomStep UtilitiesGG)

add_executable(traj GenerateTrajectories.cpp)
target_link_libraries(traj RandomStep UtilitiesGG)

#add_executable(OdeInverse OdeInverse.cpp)
#target_link_libraries(OdeInverse RandomStep UtilitiesGG MCMC)

#add_executable(OdeInverseEmb OdeInverseEmbedded.cpp)
#target_link_libraries(OdeInverseEmb RandomStep UtilitiesGG MCMC)

add_executable(HeatInverse HeatInitCond.cpp)
target_link_libraries(HeatInverse RandomStep UtilitiesGG MCMC)

add_executable(measure InvariantMeasure.cpp)
target_link_libraries(measure RandomStep UtilitiesGG)