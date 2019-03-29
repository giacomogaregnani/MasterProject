#ifndef RKSOLVER_H
#define RKSOLVER_H

#include "Structures.hpp"
#include "Tools.hpp"

class RungeKutta {

private:
    odeDef theODE;

    Butcher RKInfo;

public:
    RungeKutta() {};

    RungeKutta(odeDef ODE,
               Butcher tableau);

    VectorXd oneStep(double h,
                     VectorXd lastValue,
                     VectorXd& param);
};

#endif // RKSOLVER_H