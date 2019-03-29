#ifndef MCMC_HPP
#define MCMC_HPP

#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <memory>

#include "Posterior.hpp"
#include "Proposals.hpp"
#include "SDEPosterior.hpp"
#include "PFPosterior.hpp"
#include <memory>

using namespace Eigen;

class MCMC {
private:
    std::shared_ptr<Proposals> proposal;
    std::shared_ptr<Posterior> posterior;
    std::uniform_real_distribution<double> probGen;
    unsigned long nMCMC;
    std::vector<VectorXd> samples;
    unsigned int sizeParam;
    VectorXd initGuess;

public:
    MCMC() = default;
    ~MCMC() = default;
    MCMC(VectorXd& initGuess,
         std::shared_ptr<Proposals> proposal,
         std::shared_ptr<Posterior> posterior,
         unsigned long nMCMC);
    void eraseSample(void);
    std::vector<VectorXd>& compute(std::default_random_engine* proposalSeed,
                                   std::default_random_engine* acceptanceSeed,
                                   bool noisy = false, bool verbose = true);
};

#endif