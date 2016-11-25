add_definitions(-DDATA_PATH="${PROJECT_SOURCE_DIR}/data")
include_directories("${PROJECT_SOURCE_DIR}/Solver")

#add_executable(mlmc mlmcTest.cpp
#	       problems.hpp problems.cpp)
#target_link_libraries(mlmc Solver)

add_executable(mcmcParam mcmcParamTest.cpp
		problems.hpp problems.cpp
		mcmcTools.hpp mcmcTools.cpp
		Metropolis.cpp sMetropolis.cpp
		mlMetropolis.cpp smlMetropolis.cpp
		gaussMetropolis.cpp testMH.cpp)
target_link_libraries(mcmcParam Solver)

#add_executable(MH plainMH.cpp
#		problems.hpp problems.cpp
#		mcmcTools.hpp mcmcTools.cpp
#		Metropolis.cpp sMetropolis.cpp
#		mlMetropolis.cpp smlMetropolis.cpp
#		gaussMetropolis.cpp testMH.cpp)
#target_link_libraries(MH Solver)

add_executable(binary testBinary.cpp)

add_executable(refSol generateRefSol.cpp
               problems.hpp problems.cpp)
target_link_libraries(refSol Solver)

add_executable(testInv testInverse.cpp)
target_link_libraries(testInv Solver)

add_executable(plotter plotTrajectories.cpp
		       problems.hpp problems.cpp)
target_link_libraries(plotter Solver)

add_executable(plotGauss plotGaussTrajectories.cpp
		problems.hpp problems.cpp)
target_link_libraries(plotGauss Solver)

add_executable(jacobian jacTest.cpp
		       problems.hpp problems.cpp)

add_executable(testVar testVarMC.cpp
		       problems.hpp problems.cpp)
target_link_libraries(testVar Solver)

add_executable(likEst likelihoodEstimator.cpp
		problems.hpp problems.cpp
		mcmcTools.hpp mcmcTools.cpp
		Metropolis.cpp sMetropolis.cpp
		mlMetropolis.cpp smlMetropolis.cpp
		gaussMetropolis.cpp testMH.cpp)
target_link_libraries(likEst Solver)

#add_executable(mlmc testMLMC.cpp
#		       problems.hpp problems.cpp)
#target_link_libraries(mlmc Solver)

#add_executable(testStab testStab.cpp
#	       problems.hpp problems.cpp)
#target_link_libraries(testStab Solver)	 	
#
#add_executable(testStabErr testStabErr.cpp
#	       problems.hpp problems.cpp)
#target_link_libraries(testStabErr Solver)
#
#add_executable(smlmc smlmcTest.cpp
#	       problems.hpp problems.cpp)
#target_link_libraries(smlmc Solver)
