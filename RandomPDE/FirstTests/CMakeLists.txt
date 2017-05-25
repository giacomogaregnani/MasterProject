add_definitions(-DDATA_PATH="${PROJECT_SOURCE_DIR}/../RungeKutta/matlabPlots/ChaoticPDE/")

add_executable(Poisson1D Poisson1D.cpp)
add_executable(Poisson1D2 RandomMatrix.cpp)
add_executable(PlotPoisson PlotPoisson.cpp)
add_executable(PlotPoisson2 PlotPoissonDirOne.cpp)
add_executable(PlotPoisson3 RandomMatrixPlot.cpp)
