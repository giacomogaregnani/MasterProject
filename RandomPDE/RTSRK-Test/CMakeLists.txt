add_definitions(-DDATA_PATH="${PROJECT_SOURCE_DIR}/matlabPlots/ODETests/")

include_directories("${PROJECT_SOURCE_DIR}/RKSolver")
include_directories("${PROJECT_SOURCE_DIR}/RKRandomStep")

add_executable(Hamiltonian Hamiltonian.cpp)
target_link_libraries(Hamiltonian RandomStep)

add_executable(Hamiltonian2 Hamiltonian2.cpp)
target_link_libraries(Hamiltonian2 RandomStep)

add_executable(MeanSquareConv MeanSquareConv.cpp)
target_link_libraries(MeanSquareConv RandomStep)

add_executable(MSE MSE.cpp)
target_link_libraries(MSE RandomStep)

add_executable(Adaptivity Adaptivity.cpp)
target_link_libraries(Adaptivity RandomStep)

add_executable(GenTraj GenerateTrajectory.cpp)
target_link_libraries(GenTraj RandomStep)

add_executable(LorenzTest LorenzTestReview.cpp)
target_link_libraries(LorenzTest RandomStep)

add_executable(DetProb DistNumProb.cpp)
target_link_libraries(DetProb RandomStep)

add_executable(InvDens InvariantDensity.cpp)
target_link_libraries(InvDens RandomStep)

add_executable(EmbVsProb EmbVsProbErrEst.cpp)
target_link_libraries(EmbVsProb RandomStep)

add_executable(AdaptivityDet AdaptivityDeterministic.cpp)
target_link_libraries(AdaptivityDet RandomStep)

add_executable(OneStepVar OneStepVariance.cpp)
target_link_libraries(OneStepVar RandomStep)

add_executable(Calibration MethodCalibration.cpp)
target_link_libraries(Calibration RandomStep)