#include "RungeKuttaSolver.hpp"

RungeKutta::RungeKutta(odeDef ODE,
                       Butcher tableau):
    theODE(ODE),
    RKInfo(tableau)
{}

VectorXd RungeKutta::oneStep(double h,
                             VectorXd lastValue,
                             VectorXd& parameters)
{
    int s = RKInfo.getStages();
    VectorXd result = lastValue;

    switch (RKInfo.getType()) {

        case EXPLICIT: {

            std::vector<VectorXd> stages(s, VectorXd(theODE.size));

            // explicit methods : k_1 = f(y_0)
            stages[0] = theODE.odeFunc(lastValue, parameters);

            for (int i = 1; i < s; i++) {
                VectorXd evaluationVec = lastValue;
                for (int j = 0; j < i + 1; j++) {
                    evaluationVec = evaluationVec + h * RKInfo.getA()(i, j) * stages[j];
                }
                stages[i] = theODE.odeFunc(evaluationVec, parameters);
            }

            for (int i = 0; i < s; i++) {
                result = result + h * RKInfo.getB()(i) * stages[i];
            }
            break;
        }

        case STABEXP: {
            VectorXd kNew = lastValue, kOld = lastValue, kCurr;
            double coeff = h / static_cast<double> (s * s);
            kNew = lastValue + theODE.odeFunc(lastValue, parameters) * coeff;
            coeff = 2.0 * coeff;

            for (int i = 2; i < s + 1; i++) {
                kCurr = theODE.odeFunc(kNew, parameters) * coeff +
                        kNew * 2.0 - kOld;
                kOld = kNew;
                kNew = kCurr;
            }
            result = kCurr;
            break;
        }

        case IMPLICIT: {
            result = newtonGeneral(theODE, lastValue, parameters,
                                   h, 1000, 1e-14, RKInfo.getA(), RKInfo.getB());
            break;
        }

        case SEPAR: {
            int d = theODE.size / 2;
            if (RKInfo.getMethod() == SYMPEULER) {
                result.segment(0, d) = result.segment(0, d) + h * theODE.odeFunc(result, parameters).segment(0, d);
                result.segment(d, d) = result.segment(d, d) + h * theODE.odeFunc(result, parameters).segment(d, d);
            } else if (RKInfo.getMethod() == STORMVER) {
                result.segment(0, d) = result.segment(0, d) + 0.5 * h * theODE.odeFunc(result, parameters).segment(0, d);
                result.segment(d, d) = result.segment(d, d) + h * theODE.odeFunc(result, parameters).segment(d, d);
                result.segment(0, d) = result.segment(0, d) + 0.5 * h * theODE.odeFunc(result, parameters).segment(0, d);
            }
        }
    }

    return result;
}