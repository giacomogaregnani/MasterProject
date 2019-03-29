add_definitions(-DDATA_PATH="${PROJECT_SOURCE_DIR}/matlabPlots/CPP/")

include_directories("${PROJECT_SOURCE_DIR}/OneDimElliptic")

add_executable(DebugTest DebugTest.cpp)
target_link_libraries(DebugTest OneDimEllipticSolver)

