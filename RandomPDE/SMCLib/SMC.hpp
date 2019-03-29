#include <iostream>
#include <iomanip>
#include <memory>
#include <RandomTimeStep.hpp>


double GaussianPDF(VectorXd& x, VectorXd& m, double sigma);

class SMC {
private:

    // The number of samples
    unsigned long nParticles;

    // The update parameter
    double a;


    // Prior knowledge on the parameter
    VectorXd priorMean;
    double priorStdDev;

    // A pseudo-random Gaussian generator
    std::normal_distribution<double> gaussian;

    // Size of the parameter
    int sizeParam;

    // The initial condition of the dynamical system
    VectorXd IC;

    // A random seed
    std::default_random_engine* generator;

    // The deterministic numerical integrator
    RungeKutta detMethod;
    double h;

    // The probabilistic numerical integrator
    RungeKuttaRandomH probMethod;

public:

    SMC() {};

    SMC(unsigned long N, double a,
        VectorXd& IC, odeDef ODE, Butcher tableau,
        double h, double p,
        std::default_random_engine* generator,
        VectorXd& priorMean, double priorStdDev);

    void updateMeanAndCovariance(std::vector<VectorXd>& vec,
                                 std::vector<double>& weights,
                                 VectorXd& mean, MatrixXd& cov);

    void compute(std::vector<VectorXd>& samples,
                 std::vector<VectorXd>& thetas,
                 std::vector<double>& weights,
                 std::vector<double>& tObs,
                 std::vector<VectorXd>& yObs,
                 double obsNoise);
};