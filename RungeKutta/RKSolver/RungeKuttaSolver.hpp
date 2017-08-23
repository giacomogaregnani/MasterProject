#ifndef RKSOLVER_H
#define RKSOLVER_H

#include "Structures.hpp"
#include "Tools.hpp"

void setProblem(odeDef* ODE);

class RungeKutta {

private:
    std::vector<double> parameters;

    odeDef theODE;

    Butcher RKInfo;

public:
    RungeKutta() {};

    RungeKutta(odeDef ODE,
               std::vector<double> paramVec,
               Butcher tableau);

    VectorXd oneStep(VectorXd lastValue, double h);
};

#endif // RKSOLVER_H