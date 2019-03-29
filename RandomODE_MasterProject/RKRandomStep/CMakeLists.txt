include_directories("${PROJECT_SOURCE_DIR}/RKSolver")

add_library(RandomStep
        RandomTimeStep.hpp
        RandomTimeStep.cpp
        AdditiveNoise.hpp
        AdditiveNoise.cpp
        )
target_link_libraries(RandomStep RungeKuttaSolver)