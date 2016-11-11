#include "Solver.hpp"

double powerMethod(MatrixXd J, double tol)
{
    int N = static_cast<int>(J.rows());
    VectorXd vPower = VectorXd::Random(N);
    VectorXd tmp = J * vPower;
    double err = tol + 1;
    double minEig = 0.0;
    while (err > tol) {
        vPower = tmp / tmp.norm();
        tmp = J * vPower;
        minEig = vPower.dot(tmp) / (vPower.dot(vPower));
        err = (J * vPower - vPower * minEig).norm();
    }

    return minEig;
}