add_definitions(-DDATA_PATH="${PROJECT_SOURCE_DIR}/matlabPlots/ParameterEstimation/")

include_directories("${PROJECT_SOURCE_DIR}/NumericalSolver")
include_directories("${PROJECT_SOURCE_DIR}/ParticleFilter")
include_directories("${PROJECT_SOURCE_DIR}/MCMC")

add_executable(testParameter testParameter.cpp
               generateObservations.hpp generateObservations.cpp
               computeEstimators.hpp computeEstimators.cpp
               computeHomogeneous.hpp computeHomogeneous.cpp)
target_link_libraries(testParameter SDESolver)

add_executable(testParameterV2 testParameterV2.cpp
        generateObservations.hpp generateObservations.cpp
        computeEstimators.hpp computeEstimators.cpp
        computeHomogeneous.hpp computeHomogeneous.cpp)
target_link_libraries(testParameterV2 SDESolver)

add_executable(testParameterV2_A testParameterV2_A.cpp
        generateObservations.hpp generateObservations.cpp
        computeEstimators.hpp computeEstimators.cpp
        computeHomogeneous.hpp computeHomogeneous.cpp)
target_link_libraries(testParameterV2_A SDESolver)

add_executable(testParameterV2_S testParameterV2_S.cpp
        generateObservations.hpp generateObservations.cpp
        computeEstimators.hpp computeEstimators.cpp
        computeHomogeneous.hpp computeHomogeneous.cpp)
target_link_libraries(testParameterV2_S SDESolver)

add_executable(testPMMH testPMMH.cpp
        generateObservations.hpp generateObservations.cpp
        computeHomogeneous.hpp computeHomogeneous.cpp)
target_link_libraries(testPMMH SDESolver MCMC)

add_executable(testPMMHMod testPMMHMod.cpp
        generateObservations.hpp generateObservations.cpp
        computeHomogeneous.hpp computeHomogeneous.cpp)
target_link_libraries(testPMMHMod SDESolver MCMC)

add_executable(testHomog testHomog.cpp
        generateObservations.hpp generateObservations.cpp
        computeHomogeneous.hpp computeHomogeneous.cpp)
target_link_libraries(testHomog SDESolver)

add_executable(testModPF testModPF.cpp
        generateObservations.hpp generateObservations.cpp
        computeHomogeneous.hpp computeHomogeneous.cpp)
target_link_libraries(testModPF SDESolver MCMC)

add_executable(driftBayesian driftBayesian.cpp
        computeEstimators.hpp computeEstimators.cpp
        generateObservations.hpp generateObservations.cpp
        computeHomogeneous.hpp computeHomogeneous.cpp)
target_link_libraries(driftBayesian SDESolver MCMC)
