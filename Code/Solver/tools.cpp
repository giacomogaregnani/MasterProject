#include "Solver.hpp"

double powerMethod(MatrixXd J, double tol, int nMax)
{
    int N = static_cast<int>(J.rows());
    VectorXd vPower = VectorXd::Random(N);
    VectorXd tmp = J * vPower;
    double err = tol + 1;
    double minEig = 0.0;
    int count = 0;
    while (err > tol && count < nMax) {
        vPower = tmp / tmp.norm();
        tmp = J * vPower;
        minEig = vPower.dot(tmp) / (vPower.dot(vPower));
        err = (J * vPower - vPower * minEig).norm();
        count++;
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

VectorXd newton(odeDef ode, VectorXd uOld, std::vector<double>& param,
                double h, int nMax, double tol)
{
    double error = tol + 1;
    int iter = 0;
    VectorXd result = uOld;
    VectorXd deltaU(uOld.size());
    VectorXd F(uOld.size());
    MatrixXd JF;
    MatrixXd I = MatrixXd::Identity(uOld.size(), uOld.size());
    while (error > tol && iter++ <= nMax) {
        JF = I - h * ode.odeJac(result, param);
        F = result - h * ode.odeFunc(result, param) - uOld;
        deltaU = -1.0 * JF.partialPivLu().solve(F);
        result = result + deltaU;
        error = deltaU.norm();
    }
    // printf("Newton converged to %e in %d iterations\n", error, iter);
    return result;
}


VectorXd newtonGeneral(odeDef ode, VectorXd uOld, std::vector<double>& param,
                       double h, int nMax, double tol, MatrixXd A, VectorXd b)
{
    long int s = b.size(), d = uOld.size();
    double error = tol + 1;
    int iter = 0;

    VectorXd result = uOld;
    VectorXd K = uOld.replicate(s, 1);
    MatrixXd I = MatrixXd::Identity(s * d, s * d);
    MatrixXd JF(s*d, s*d), evaluatedJacSmallF(d, d);
    VectorXd F(s*d), deltaK(s * d), RHS(s*d);

    evaluatedJacSmallF = ode.odeJac(uOld, param);
    JF = I - kroneckerProduct(A, evaluatedJacSmallF) * h;
    auto JFfact = JF.fullPivLu();

    while (error > tol && iter++ <= nMax) {
        // Assemble the RHS
        for (int i = 0; i < s; i++) {
            F.segment(d * i, d) = VectorXd::Zero(d);
            for (int j = 0; j < s; j++) {
                F.segment(d * i, d) += h * A(i, j) * ode.odeFunc(uOld + K.segment(d * j, d), param);
            }
        }
        RHS = F - K;

        // Newton step
        deltaK = JFfact.solve(RHS);
        K = K + deltaK;
        error = deltaK.norm();
    }
    // printf("Newton converged to %e in %d iterations\n", error, iter);


    // Assemble the result
    for (int i = 0; i < s; i++) {
        result += h * b(i) * ode.odeFunc(uOld + K.segment(d * i, d), param);
    }

    return result;
}

// BUTCHER TABLES TOOLS
Butcher::Butcher(implicitMethods chosenMethod, int s)
{
    method = chosenMethod;
    stages = s;
    A.resize(2, 2);
    b.resize(2);

    switch (method) {
        case GAUSS:
            A(0, 0) = 0.25;
            A(0, 1) = 0.25 - sqrt(3.0) / 6.0;
            A(1, 0) = 0.25 + sqrt(3.0) / 6.0;
            A(1, 1) = 0.25;
            b(0) = 0.5;
            b(1) = 0.5;
            break;
        case LOBATTO:
            A(0, 0) = 0.5;
            A(0, 1) = 0.0;
            A(1, 0) = 0.5;
            A(1, 1) = 0.0;
            b(0) = 0.5;
            b(1) = 0.5;
            break;
        case RADAU:
            A(0, 0) = 0.25;
            A(0, 1) = -0.25;
            A(1, 0)  = 0.25;
            A(1, 1) = 5.0 / 12.0;
            b(0) = 0.25;
            b(1) = 0.75;
            break;
    }
}

MatrixXd Butcher::getA(void)
{
    return A;
}

VectorXd Butcher::getB(void)
{
    return b;
}

int Butcher::getStages(void)
{
    return stages;
}

