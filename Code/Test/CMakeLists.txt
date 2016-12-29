add_definitions(-DDATA_PATH="${PROJECT_SOURCE_DIR}/data")
include_directories("${PROJECT_SOURCE_DIR}/Solver")

add_executable(testMerge testMerge.cpp)

add_executable(LorenzTraj LorenzPlotTrajectories.cpp
		problems.cpp problems.hpp)
target_link_libraries(LorenzTraj Solver)

add_executable(nonZeroVar nonZeroVar.cpp
		problems.cpp problems.hpp)
target_link_libraries(nonZeroVar Solver)

#add_executable(MH plainMH.cpp
#		problems.hpp problems.cpp
#		mcmcTools.hpp mcmcTools.cpp
#		Metropolis.cpp sMetropolis.cpp
#		mlMetropolis.cpp smlMetropolis.cpp
#		gaussMetropolis.cpp testMH.cpp)
#target_link_libraries(MH Solver)

add_executable(testBE testBackEul.cpp
		       problems.hpp problems.cpp)
target_link_libraries(testBE Solver)

add_executable(binary testBinary.cpp)

add_executable(refSol generateRefSol.cpp
               problems.hpp problems.cpp)
target_link_libraries(refSol Solver)

add_executable(plotter plotTrajectories.cpp
		       problems.hpp problems.cpp)
target_link_libraries(plotter Solver)

add_executable(plotGauss plotGaussTrajectories.cpp
		problems.hpp problems.cpp)
target_link_libraries(plotGauss Solver)

add_executable(testVar testVarMC.cpp
		       problems.hpp problems.cpp)
target_link_libraries(testVar Solver)

add_executable(plainMH plainMH.cpp
		       mcmcTools.hpp mcmcTools.cpp
		       testMH.cpp)

add_executable(plainMHconv plainMHConv.cpp
		mcmcTools.hpp mcmcTools.cpp
		testMHconv.cpp)

add_executable(invMeas invMeasure.cpp
		problems.hpp problems.cpp)
target_link_libraries(invMeas Solver)

#add_executable(likEst likelihoodEstimator.cpp
#		problems.hpp problems.cpp
#		mcmcTools.hpp mcmcTools.cpp
#		Metropolis.cpp sMetropolis.cpp
#       MetropolisConv.cpp
#		#mlMetropolis.cpp
#		smlMetropolis.cpp
#		MetropolisDet.cpp
#		gaussMetropolis.cpp testMH.cpp)
#target_link_libraries(likEst Solver)

add_executable(mcmcConv testMCMCconvergence.cpp
		problems.hpp problems.cpp
		mcmcTools.hpp mcmcTools.cpp
		Metropolis.cpp sMetropolis.cpp
		MetropolisConv.cpp
		#mlMetropolis.cpp
		impMetropolisConv.cpp
		smlMetropolis.cpp
		MetropolisDet.cpp
		gaussMetropolis.cpp testMH.cpp)
target_link_libraries(mcmcConv Solver)

add_executable(mcmcConvFix MCMCconvFix.cpp
		problems.hpp problems.cpp
		mcmcTools.hpp mcmcTools.cpp
		Metropolis.cpp sMetropolis.cpp
		MetropolisConv.cpp
		#mlMetropolis.cpp
		impMetropolisConv.cpp
		smlMetropolis.cpp
		MetropolisDet.cpp
		gaussMetropolis.cpp testMH.cpp)
target_link_libraries(mcmcConvFix Solver)

add_executable(sigmaEst sigmaEst.cpp
		problems.hpp problems.cpp
		mcmcTools.hpp mcmcTools.cpp
		Metropolis.cpp sMetropolis.cpp
		MetropolisConv.cpp
		#mlMetropolis.cpp
		impMetropolisConv.cpp
		smlMetropolis.cpp
		MetropolisDet.cpp
		gaussMetropolis.cpp testMH.cpp)
target_link_libraries(sigmaEst Solver)

add_executable(weakConv weakConvergence.cpp
		problems.hpp problems.cpp)
target_link_libraries(weakConv Solver)