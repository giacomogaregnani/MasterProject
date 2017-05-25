include_directories("${PROJECT_SOURCE_DIR}/MCMCLib")
add_definitions(-DDATA_PATH="${PROJECT_SOURCE_DIR}/matlabPlots/")

add_executable(MCMCTest testMCMC.cpp)
target_link_libraries(MCMCTest MCMC)