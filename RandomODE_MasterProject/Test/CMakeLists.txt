add_definitions(-DDATA_PATH="${PROJECT_SOURCE_DIR}/matlabScripts/")

include_directories("${PROJECT_SOURCE_DIR}/RKSolver")
include_directories("${PROJECT_SOURCE_DIR}/RKRandomStep")

add_executable(example GenerateTrajectories.cpp)
target_link_libraries(example RandomStep)

add_executable(strongConv StrongConvergence.cpp)
target_link_libraries(strongConv RandomStep)

add_executable(weakConv WeakConvergence.cpp)
target_link_libraries(weakConv RandomStep)