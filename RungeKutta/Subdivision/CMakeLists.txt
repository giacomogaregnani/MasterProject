include_directories("${PROJECT_SOURCE_DIR}/RKSolver")

add_library(Subdivision
        Subdivision.hpp
        Subdivision.cpp
        Box.hpp
        Box.cpp
        Basics.hpp
        )

target_link_libraries(Subdivision RungeKuttaSolver)