include_directories("${PROJECT_SOURCE_DIR}/RKRandomStep")
include_directories("${PROJECT_SOURCE_DIR}/RKSolver")

add_library(SMC
        SMC.hpp
        SMC.cpp
        SMCDet.hpp
        SMCDet.cpp
        )

target_link_libraries(SMC RandomStep)