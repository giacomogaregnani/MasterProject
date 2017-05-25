#include "MCMC.hpp"


MCMC::MCMC(std::normal_distribution<double>::param_type& proposalParam,
           ProbDensFunc *prior, LikelihoodFunc *likelihood,
           std::vector<VectorXd>& observations,
           unsigned long nMCMC, bool RAM, VectorXd initGuess,
           double desiredAlpha):
        prior(prior),
        likelihood(likelihood),
        nMCMC(nMCMC),
        RAM(RAM),
        observations(observations)
{
    probGen.param(std::uniform_real_distribution<double>::param_type(0.0, 1.0));
    samples = {initGuess};
    sizeParam = initGuess.size();
    if (RAM) {
        RAMUpdate.init(proposalParam.stddev(), desiredAlpha, initGuess.size());
        proposal.param(std::normal_distribution<double>::param_type(0.0, 1.0));
    } else {
        proposal.param(proposalParam);
    }
}

std::vector<VectorXd>& MCMC::compute(std::default_random_engine* generator)
{
    double newPosterior, newLikelihood, newPrior,
           oldPosterior, oldLikelihood, oldPrior, alpha, currBest;
    unsigned long accRatio = 0;

    VectorXd increment(sizeParam), proposedValue(sizeParam);

    for (unsigned long i = 0; i < nMCMC; i++) {

        // generate new guess
        for (unsigned int j = 0; j < sizeParam; j++) {
            increment(j) = proposal(*generator);
        }
        if (RAM) {
            proposedValue = samples.back() + RAMUpdate.getS() * increment;
        } else {
            proposedValue = samples.back() + increment;
        }

        // evaluate posterior on new guess
        newLikelihood = likelihood->LIK(proposedValue, observations);
        newPrior = prior->PDF(proposedValue);
        newPosterior = newPrior + newLikelihood;

        // re-evaluate posterior on old guess
        oldLikelihood = likelihood->LIK(samples.back(), observations);
        oldPrior = prior->PDF(samples.back());
        oldPosterior = oldPrior + oldLikelihood;

        // compute probability
        alpha = newPosterior - oldPosterior;
        alpha = std::min(1.0, exp(alpha));

        // Update chain
        if (probGen(*generator) < alpha) {
            samples.push_back(proposedValue);
            oldPosterior = newPosterior;
            ++accRatio;
            currBest = newLikelihood;
        } else {
            samples.push_back(samples.back());
        }

        if (RAM) {
            RAMUpdate.update(increment, alpha);
        }

        if (i % 100 == 0) {
            std::cout << "Completed " << i << " iterations out of " << nMCMC << std::endl
                      << "Parameter value = " << samples.back().transpose() << std::endl
                      << "log-likelihood = " <<  currBest << std::endl;
        }
    }

    std::cout << "=====================" << std::endl
              << "MCMC END WITH "
              << std::fixed << std::setprecision(3)
              << static_cast<double>(accRatio) / nMCMC
              << " ACCEPTANCE RATIO." << std::endl
              << "=====================" << std::endl;

    return samples;
}

