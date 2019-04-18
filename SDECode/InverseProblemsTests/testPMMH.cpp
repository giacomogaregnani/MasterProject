#include <fstream>
#include <iostream>
#include <iomanip>
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

    double T = 20;
    unsigned int N = 4000;

    VectorXd tmpParam(3);
    tmpParam(0) = 0.1;  // Epsilon
    tmpParam(1) = 1.0;  // True multiscale alpha
    tmpParam(2) = 0.5;  // True multiscale betainv
    // tmpParam(3) = std::log(3*T / N); // log of the sampling period

    std::ofstream output(DATA_PATH + std::string("testHomo2.txt"), std::ofstream::out | std::ofstream::trunc);
    std::ofstream outputSol(DATA_PATH + std::string("testHomoSol2.txt"), std::ofstream::out | std::ofstream::trunc);
    std::ofstream outputFilter(DATA_PATH + std::string("testHomoFil2.txt"), std::ofstream::out | std::ofstream::trunc);

    // Compute coefficients of the homogenised equation
    std::vector<double> homCoeffs = computeHomCoeffs(tmpParam, (2.0*M_PI), &V1);
    output << homCoeffs[0] << "\t" << homCoeffs[1] << std::endl;
    homoDiff = homCoeffs[1];

    // ============= Transform for positiveness ============= //
    VectorXd param = tmpParam;
    param(1) = std::log(param(1));
    param(2) = std::log(param(2));
    // ====================================================== //

    // Initialize structures for the inverse problem
    unsigned long M = 100, nMCMC = 1001;
    double noise = 1e-2;
    double IC = 0.0;
    std::random_device dev;
    std::default_random_engine noiseSeed{1};
    //std::default_random_engine noiseSeed{dev()};
    std::default_random_engine proposalSeed{(unsigned int) time(nullptr)+1};
    std::default_random_engine acceptanceSeed{(unsigned int) time(nullptr)+2};
    std::normal_distribution<double> noiseDistribution(0.0, noise);
    bool IS = true;

    // Generate and perturb observations
    auto x = generateObservations1D(sde, IC, param, T, N);
    for (auto const &itSol : x) {
        outputSol << std::fixed << std::setprecision(5) << itSol << "\t";
    }
    outputSol << std::endl;
    outputSol << x[0] << "\t";
    for (unsigned long i = 1; i < x.size(); i++) {
        x[i] += noiseDistribution(noiseSeed);
        outputSol << x[i] << "\t";
    }
    outputSol << std::endl;
    // Compute a homogeneous path with the same Brownian for comparison
    VectorXd homoParam(3);
    homoParam(0) = 0.1;
    homoParam(1) = std::log(homCoeffs[0]);
    homoParam(2) = std::log(homCoeffs[1]);
    auto xHom = generateObservations1D(sdeHomo, IC, homoParam, T, N);
    for (auto const &itSol : xHom) {
        outputSol << std::fixed << std::setprecision(5) << itSol << "\t";
    }
    outputSol << std::endl;

    // Initial parameter guess
    VectorXd initGuess = VectorXd::Zero(param.size());
    // initGuess(3) = param(3);

    // Signal ratio
    std::vector<unsigned int> ratioVec = {1};

    // Inverse problem
    for (auto const &ratio : ratioVec) {
        std::shared_ptr<Posterior> posterior;
        posterior = std::make_shared<PFPosterior>(x, T, IC, ratio, noise, sdeHomo, param(0), M, IS);
        std::shared_ptr<Proposals> proposal;
        proposal = std::make_shared<Proposals>(1e-1);
        MCMC mcmc(initGuess, proposal, posterior, nMCMC);
        auto sample = mcmc.compute(&proposalSeed, &acceptanceSeed);
        for (auto const &itSample : sample)
            output << itSample.transpose() << std::endl;

        // Compute the mean and write on file the pushed mean
        VectorXd meanPosterior = VectorXd::Zero(3);
        for (auto const &it : sample)
            meanPosterior += it;
        meanPosterior /= sample.size();

        /*
        ParFil particleFilter(x, T, IC, ratio, noise, sdeHomo, param(0), M);
        particleFilter.computeDiffBridge(meanPosterior);
        auto XPost = particleFilter.getX();
        for (auto const &itSol : XPost) {
            for (auto const &itit : itSol)
                outputFilter << itit << "\t";
            outputFilter << std::endl;
        }
        */

        M = static_cast<unsigned long>(M * 1.2);
    }

    output.close();
    outputSol.close();
    outputFilter.close();
    return 0;
}