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

MatrixXd triInv(MatrixXd A, int size)
{
    int i, j, k;
    double tmp;

    for (j = 0; j < size; j++) {
        A(j, j) = 1.0 / A(j, j);
        for (i = j + 1; i < size; i++) {
            tmp = 0.0;
            for (k = j; k < i; k++) {
                tmp -= A(i, k) * A(k, j);
            }
            A(i, j) = tmp / A(i, i);
        }
    }

    return A;
}
