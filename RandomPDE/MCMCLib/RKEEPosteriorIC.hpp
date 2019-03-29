#ifndef RKEEPOSTERIORIC_HPP
#define RKEEPOSTERIORIC_HPP

#include <Eigen/Dense>
#include <RungeKuttaSolver.hpp>
#include "Proposals.hpp"
#include "Posterior.hpp"
#include "KarhunenLoeve.hpp"

using namespace Eigen;

class RKEEPosteriorIC : public Posterior {
private:

    odeDef ODE;

    RungeKutta RKSolver;

    double h;

    std::vector<VectorXd> observations;

    std::vector<double> tObs;

    double noise;

    MatrixXd noiseMatrix;

    LLT<MatrixXd> noiseMatrixLLT;

    bool isIdentity;

    std::vector<VectorXd> modErrors;

public:

    RKEEPosteriorIC() {};

    RKEEPosteriorIC(double h,
                    std::vector<VectorXd>& observations,
                    std::vector<VectorXd>& modErrors,
                    std::vector<double>& tObs,
                    double noise,
                    odeDef ODE,
                    Butcher tableau);

    RKEEPosteriorIC(double h,
                    std::vector<VectorXd>& observations,
                    std::vector<VectorXd>& modErrors,
                    std::vector<double>& tObs,
                    MatrixXd& noise,
                    odeDef ODE,
                    Butcher tableau);

    double computePosterior(VectorXd& theta);

    VectorXd coeffToField(VectorXd& theta) {};
};


#endif