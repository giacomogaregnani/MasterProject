add_definitions(-DDATA_PATH_SRC="${PROJECT_SOURCE_DIR}/data")

add_library(Solver Solver.hpp
		   Structures.hpp
		   ProbMethod.cpp
		   sProbMethod.cpp
		   impProbMethod.cpp
		   smlmc.cpp
		   RungeKutta.hpp RungeKutta.cpp
		   sRungeKutta.hpp sRungeKutta.cpp
		   ThirdOrderGauss.cpp tools.cpp)