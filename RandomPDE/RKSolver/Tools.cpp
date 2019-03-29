#include "Tools.hpp"

VectorXd newtonGeneral(odeDef ode, VectorXd uOld, VectorXd& param,
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

    // Assemble the result
    for (int i = 0; i < s; i++) {
        result += h * b(i) * ode.odeFunc(uOld + K.segment(d * i, d), param);
    }

    return result;
}