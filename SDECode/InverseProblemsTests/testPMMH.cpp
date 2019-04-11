#include <fstream>
#include <iostream>
#include "generateObservations.hpp"
#include "computeHomogeneous.hpp"
#include <MCMC.hpp>

// Remark: epsilon = p(0)

double homoDiff;

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
    return -1.0 * (std::exp(p(1)) * gradV0(x) + 1 / p(0) * gradV1(x / p(0)));
}

double homoDrift(double x, VectorXd& p)
{
    return -1.0 * std::exp(p(1)) * gradV0(x);
}

double diffusion(double x, VectorXd &p)
{
    // return std::sqrt(2.0 * 0.5);
    return std::sqrt(2.0 * std::exp(p(2)));
}

double homoDiffusion(double x, VectorXd& p)
{
    return std::sqrt(2.0 * homoDiff);
}


int main(int argc, char* argv[])
{
    oneDimSde sde;
    sde.drift = &multiDrift;
    sde.diffusion = &diffusion;
    oneDimSde sdeHomo;
    sdeHomo.drift = &homoDrift;
    sdeHomo.diffusion = &diffusion;

    VectorXd tmpParam(3);
    tmpParam(0) = 0.1;  // Epsilon
    tmpParam(1) = 1.0;  // True multiscale alpha
    tmpParam(2) = 0.5;  // True multiscale betainv

    double T = 40;
    unsigned int N = 4000;

    std::ofstream output(DATA_PATH + std::string("testHomo.txt"), std::ofstream::out | std::ofstream::trunc);
    std::ofstream outputSol(DATA_PATH + std::string("testHomoSol.txt"), std::ofstream::out | std::ofstream::trunc);

    // Compute coefficients of the homogenised equation
    std::vector<double> homCoeffs = computeHomCoeffs(tmpParam, (2.0*M_PIf64), &V1);
    output << homCoeffs[0] << "\t" << homCoeffs[1] << std::endl;
    homoDiff = homCoeffs[1];

    // ============= Transform for positiveness ============= //
    VectorXd param = tmpParam;
    param(1) = std::log(param(1));
    param(2) = std::log(param(2));
    // ====================================================== //

    // Initialize structures for the inverse problem
    unsigned long M = 200, nMCMC = 5000;
    double noise = 1e-2;
    double IC = 0.0;
    VectorXd initGuess = VectorXd::Zero(param.size());
    initGuess(0) = param(0);
    std::default_random_engine noiseSeed{1};
    std::default_random_engine proposalSeed{(unsigned int) time(nullptr)+1};
    std::default_random_engine acceptanceSeed{(unsigned int) time(nullptr)+2};
    std::normal_distribution<double> noiseDistribution(0.0, noise);

    // Signal ratio
    std::vector<unsigned int> ratioVec = {200};

    // Generate and perturb observations
    auto x = generateObservations1D(sde, IC, param, T, N);
    for (auto const &itSol : x) {
        outputSol << itSol << "\t";
    }
    outputSol << std::endl;
    outputSol << x[0] << "\t";
    for (unsigned long i = 1; i < x.size(); i++) {
        x[i] += noiseDistribution(noiseSeed);
        outputSol << x[i] << "\t";
    }
    outputSol << std::endl;

    // Inverse problem
    for (auto const &ratio : ratioVec) {
        std::shared_ptr<Posterior> posterior;
        posterior = std::make_shared<PFPosterior>(x, T, IC, ratio, noise, sdeHomo, param(0), M, true);
        std::shared_ptr<Proposals> proposal;
        proposal = std::make_shared<Proposals>(1e-1);
        MCMC mcmc(initGuess, proposal, posterior, nMCMC);
        auto sample = mcmc.compute(&proposalSeed, &acceptanceSeed);
        for (auto const &itSample : sample)
            output << itSample.transpose() << std::endl;
    }
    output.close();

    return 0;
}