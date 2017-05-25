#ifndef TOOLS_H
#define TOOLS_H

#include "Structures.hpp"
#include <unsupported/Eigen/KroneckerProduct>
#include <iostream>

class Poisson {
private:
    static MatrixXd A;

public:
    static VectorXd initialCond;
    static VectorXd poissonTest(VectorXd argument, std::vector<double>& param);
    static MatrixXd poissonTestJ(VectorXd argument, std::vector<double>& param);
    static VectorXd poissonExact(std::vector<double>& param, double t);
    Poisson(int N);
};


VectorXd newtonGeneral(odeDef ode, VectorXd uOld, std::vector<double>& param,
                       double h, int nMax, double tol, MatrixXd A, VectorXd b);


#endif
