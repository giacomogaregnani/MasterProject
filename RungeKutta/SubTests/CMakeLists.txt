add_definitions(-DDATA_PATH="${PROJECT_SOURCE_DIR}/matlabPlots/")

include_directories("${PROJECT_SOURCE_DIR}/RKSolver")
include_directories("${PROJECT_SOURCE_DIR}/Subdivision")

add_executable(Basic BasicTest.cpp)
target_link_libraries(Basic Subdivision)

