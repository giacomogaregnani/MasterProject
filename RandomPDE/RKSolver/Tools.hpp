#ifndef TOOLS_H
#define TOOLS_H

#include "Structures.hpp"
#include <unsupported/Eigen/KroneckerProduct>
#include <iostream>

class Heat {
private:
    static SparseMatrix<double> A;

public:
    static VectorXd initialCond;
    static VectorXd heatF(VectorXd argument, std::vector<double> &param);
    static MatrixXd heatJ(VectorXd argument, std::vector<double> &param);
    static VectorXd poissonExact(std::vector<double>& param, double t);
    Heat(int N);
};


VectorXd newtonGeneral(odeDef ode, VectorXd uOld, VectorXd& param,
                       double h, int nMax, double tol, MatrixXd A, VectorXd b);


#endif
