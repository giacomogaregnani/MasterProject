include_directories("${PROJECT_SOURCE_DIR}/NumericalSolver")

add_library(ParFil
        ParFil.hpp
        ParFil.cpp
        )

target_link_libraries(ParFil SDESolver)