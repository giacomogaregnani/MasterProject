#ifndef RKSOLVER_H
#define RKSOLVER_H

#include "Structures.hpp"
#include "Tools.hpp"


class RungeKutta {

private:
    LinearEquation theODE;

    Butcher RKInfo;

    FullPivLU<MatrixXd> ImplicitMatrix;

    MatrixXd RHSMatrix;

public:
    RungeKutta() {};

    RungeKutta(LinearEquation ODE,
               Butcher tableau);

    VectorXd oneStep(double h,
                     VectorXd lastValue);
};

#endif // RKSOLVER_H