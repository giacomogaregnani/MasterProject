include_directories("${PROJECT_SOURCE_DIR}/NumericalSolver")

add_library(ParFil
        ParFilLib.hpp
        ParFil.hpp
        ParFil.cpp
        ParFilMod.hpp
        ParFilMod.cpp
        ParFilTools.hpp
        ParFilTools.cpp
        )

target_link_libraries(ParFil SDESolver)