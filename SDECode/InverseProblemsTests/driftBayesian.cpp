#include <fstream>
#include <iostream>
#include <iomanip>
#include <MCMC.hpp>
#include "generateObservations.hpp"
#include "computeHomogeneous.hpp"
#include "computeEstimators.hpp"

// Remark: epsilon = p(0)

double gradV0(double x)
{
    return x;
}

double V1(double x)
{
    return std::cos(x);
}

double gradV1(double x)
{
    return -std::sin(x);
}

double multiDrift(double x, VectorXd& p)
{
    return -1.0 * (p(1) * gradV0(x) + 1 / p(0) * gradV1(x / p(0)));
}

double diffusion(double x, VectorXd &p)
{
    return std::sqrt(2.0 * p(2));
}

int main(int argc, char* argv[])
{
    oneDimSde sde{&multiDrift, &diffusion};
    double IC = 0.0;

    VectorXd param(3);
    param(0) = 0.1;
    param(1) = 1.0; // True multiscale alpha
    param(2) = 0.5; // True multiscale sigma

    std::ofstream output(DATA_PATH + std::string("driftBayesian.txt"), std::ofstream::out | std::ofstream::trunc);

    // Compute coefficients of the homogenised equation
    std::vector<double> homCoeffs = computeHomCoeffs(param, (2.0*M_PI), &V1);
    output << homCoeffs[0] << std::endl;

    // Compute the estimators for different value of epsilon
    std::vector<double> avg;

    // Refer to the Caltech notes
    double beta = 2.0;
    double zeta = 2.0;
    double gamma = 5.0;

    std::random_device dev;
    std::default_random_engine proposalSeed(dev());
    std::default_random_engine acceptanceSeed(dev());

    double eps = param(0);

    auto h = std::pow(eps, beta);
    auto N = static_cast<unsigned int>(std::round(std::pow(eps, -gamma)));
    auto T = h * N;
    auto delta = static_cast<unsigned int>(std::round(std::pow(eps, -zeta)));

    auto trajectory = generateObservations1D(sde, IC, param, T, N, dev());
    avg = averageSequence(trajectory, delta);

    auto obs = avg;

    // First Compute MLE
    double AHat = estimateA(obs, h, &gradV0);
    output << AHat << std::endl;

    // Then sample from posterior
    std::shared_ptr<Posterior> posterior;
    posterior = std::make_shared<NLPosterior>(obs, h, &gradV0, param);
    std::shared_ptr<Proposals> proposal;
    proposal = std::make_shared<Proposals>(5e-2);

    VectorXd initGuess = VectorXd::Zero(1);

    MCMC mcmc(initGuess, proposal, posterior, 2000);
    auto sample = mcmc.compute(&proposalSeed, &acceptanceSeed);

    // Sample from homogenized posterior
    // TODO

    for (auto it : sample)
        output << it << std::endl;

    output.close();
    return 0;
}