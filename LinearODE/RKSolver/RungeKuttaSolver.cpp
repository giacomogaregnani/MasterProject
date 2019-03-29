#include "RungeKuttaSolver.hpp"
#include <unsupported/Eigen/KroneckerProduct>

RungeKutta::RungeKutta(LinearEquation ODE,
                       Butcher tableau):
    theODE(ODE),
    RKInfo(tableau)
{
    if (RKInfo.getType() == IMPLICIT) {
        MatrixXd Is = MatrixXd::Identity(RKInfo.getStages(), RKInfo.getStages());
        MatrixXd Id = MatrixXd::Identity(theODE.getSize(), theODE.getSize());
        ImplicitMatrix = FullPivLU<MatrixXd>(kroneckerProduct(theODE.getM(), Is) -
                                             kroneckerProduct(theODE.getA(), RKInfo.getA()));
        RHSMatrix = kroneckerProduct(theODE.getM(), Is);
    }

}

VectorXd RungeKutta::oneStep(double h, VectorXd lastValue)
{
    int s = RKInfo.getStages();
    VectorXd result = lastValue;

    switch (RKInfo.getType()) {

        case EXPLICIT: {
            std::vector<VectorXd> stages(s, VectorXd(theODE.getSize()));

            // explicit methods : k_1 = f(y_0)
            stages[0] = theODE.f(lastValue);

            for (int i = 1; i < s; i++) {
                VectorXd evaluationVec = lastValue;
                for (int j = 0; j < i + 1; j++) {
                    evaluationVec = evaluationVec + h * RKInfo.getA()(i, j) * stages[j];
                }
                stages[i] = theODE.f(evaluationVec);
            }

            for (int i = 0; i < s; i++) {
                result = result + h * RKInfo.getB()(i) * stages[i];
            }
            break;
        }

        case STABEXP: {
            VectorXd kNew = lastValue, kOld = lastValue, kCurr;
            double coeff = h / static_cast<double> (s * s);
            kNew = lastValue + theODE.f(lastValue) * coeff;
            coeff = 2.0 * coeff;

            for (int i = 2; i < s + 1; i++) {
                kCurr = theODE.f(kNew) * coeff +
                        kNew * 2.0 - kOld;
                kOld = kNew;
                kNew = kCurr;
            }
            result = kCurr;
            break;
        }

        case IMPLICIT: {
            result = newtonGeneral(theODE, lastValue, h, 1000, 1e-14, RKInfo.getA(), RKInfo.getB());
            break;
        }
    }

    return result;
}