#include <RungeKuttaSolver.hpp>

struct forwardMapStruct {
    forwardMapStruct() {}

    odeDef ODE;
    std::vector<double> paramVec;
    Butcher tableau;
    double h;
};

