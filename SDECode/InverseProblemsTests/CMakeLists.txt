add_definitions(-DDATA_PATH="${PROJECT_SOURCE_DIR}/matlabPlots/ParameterEstimation/")

include_directories("${PROJECT_SOURCE_DIR}/NumericalSolver")
include_directories("${PROJECT_SOURCE_DIR}/ParticleFilter")
include_directories("${PROJECT_SOURCE_DIR}/MCMC")

add_executable(testParameter testParameter.cpp
               generateObservations.hpp generateObservations.cpp
               computeEstimators.hpp computeEstimators.cpp
               computeHomogeneous.hpp computeHomogeneous.cpp)
target_link_libraries(testParameter SDESolver)

add_executable(testPMMH testPMMH.cpp
        generateObservations.hpp generateObservations.cpp
        computeHomogeneous.hpp computeHomogeneous.cpp)
target_link_libraries(testPMMH SDESolver MCMC)

add_executable(testModErr testModErr.cpp
        generateObservations.hpp generateObservations.cpp)
target_link_libraries(testModErr SDESolver MCMC)

add_executable(testPMMHMod testPMMHModErr.cpp
        generateObservations.hpp generateObservations.cpp
        computeHomogeneous.hpp computeHomogeneous.cpp)
target_link_libraries(testPMMHMod SDESolver MCMC)

add_executable(testPMMHHom testPMMHHom.cpp
        generateObservations.hpp generateObservations.cpp
        computeHomogeneous.hpp computeHomogeneous.cpp)
target_link_libraries(testPMMHHom SDESolver MCMC)

add_executable(testHomog TestHomog.cpp
        generateObservations.hpp generateObservations.cpp
        computeHomogeneous.hpp computeHomogeneous.cpp)
target_link_libraries(testHomog SDESolver)