#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <memory>

#include "FEMPosterior.hpp"
#include "FEMPosteriorFin.hpp"
#include "FEMProbPosterior.hpp"
#include "FEMProbPosteriorMulti.hpp"
#include "FEMProbPosteriorFin.hpp"
#include "RKPosterior.hpp"
#include "RKProbPosterior.hpp"
#include "RKProbPosteriorFixed.hpp"
#include "RKPosteriorIC.hpp"
#include "RKProbPosteriorIC.hpp"
#include "RKEEPosteriorIC.hpp"
#include "RAM.hpp"
#include "Posterior.hpp"
#include "RKAddNoisePosteriorIC.hpp"

using namespace Eigen;

class MCMC {
private:

    // The proposal distribution
    Proposals* proposal;

    // The posterior distribution
    Posterior* posterior;

    // The probability generator
    std::uniform_real_distribution<double> probGen;

    // The number of samples
    unsigned long nMCMC;

    // The generated samples
    std::vector<VectorXd> samples;
    unsigned int sizeParam;

    // The initial guess
    VectorXd initGuess;

public:

    MCMC() = default;

    MCMC(VectorXd& initGuess,
         Proposals* proposal,
         Posterior* posterior,
         unsigned long nMCMC);

    void eraseSample(void);

    std::vector<VectorXd>& compute(std::default_random_engine* proposalSeed,
                                   std::default_random_engine* acceptanceSeed,
                                   bool noisy = false, bool verbose = true);
};