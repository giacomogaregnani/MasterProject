#include "MCMC.hpp"

MCMC::MCMC(VectorXd& initGuess,
           std::normal_distribution<double>::param_type& proposalParam,
           double (*post) (VectorXd& theta),
           unsigned long nMCMC,
           bool RAM, double desiredAlpha):
        posterior(post),
        nMCMC(nMCMC),
        RAM(RAM)
{
    probGen.param(std::uniform_real_distribution<double>::param_type(0.0, 1.0));
    samples = {initGuess};
    sizeParam = static_cast<unsigned int> (initGuess.size());
    if (RAM) {
        RAMUpdate.init(proposalParam.stddev(), desiredAlpha, sizeParam);
        proposal.param(std::normal_distribution<double>::param_type(0.0, 1.0));
    } else {
        proposal.param(proposalParam);
    }
}

std::vector<VectorXd>& MCMC::compute(std::default_random_engine* generator, bool noisy)
{
    double newPosterior, oldPosterior, alpha;
    VectorXd currEstimate;
    currEstimate = samples.back();
    unsigned long accRatio = 0;

    VectorXd increment(sizeParam), proposedValue(sizeParam);

    // Compute posterior on first guess
    oldPosterior = posterior(samples.back());

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
        newPosterior = posterior(proposedValue);

        // For the noisy algorithm, re-evaluate on old guess
        if (noisy) {
            oldPosterior = posterior(samples.back());
        }

        // compute probability
        alpha = newPosterior - oldPosterior;
        alpha = std::min(0.0, alpha);

        // Update chain
        if (std::log(probGen(*generator)) < alpha) {
            samples.push_back(proposedValue);
            oldPosterior = newPosterior;
            ++accRatio;
        } else {
            samples.push_back(samples.back());
        }

        // Update current estimate
        currEstimate = currEstimate * static_cast<double>(i + 1) / (i + 2)
                       + samples.back() * 1.0 / (i + 2);

        if (RAM) {
            RAMUpdate.update(increment, exp(alpha));
        }

        if (i % 1000 == 0) {
            std::cout << "Completed " << i << " iterations out of " << nMCMC << std::endl
                      << "Current estimate = " << currEstimate.transpose() << std::endl
                      << "Last accepted value = " << samples.back().transpose() << std::endl
                      << "Current acc. ratio = " << 1.0 / (i + 1) * accRatio << std::endl
                      << "Current posterior = " << oldPosterior << std::endl
                      << "==============================" << std::endl;
        }
    }

    // Take 10% burn-in out
    unsigned long burnIn = nMCMC / 10;
    samples.erase(samples.begin(), samples.begin() + burnIn);

    // Compute MCMC estimate
    VectorXd finalEstimate = VectorXd::Zero(currEstimate.size());
    for (auto it : samples) {
        finalEstimate += it;
    }
    finalEstimate /= (nMCMC - burnIn);

    std::cout << std::fixed << std::setprecision(3)
              << "=====================" << std::endl
              << "MCMC ESTIMATE:" << std::endl
              << finalEstimate.transpose() << std::endl
              << " ACCEPTANCE RATIO: " << static_cast<double>(accRatio) / nMCMC << std::endl;

    if (RAM) {
        std::cout << "FINAL COVARIANCE MATRIX (proposal)" << std::endl
                  << RAMUpdate.getS() * RAMUpdate.getS().transpose() << std::endl
                  << "=====================" << std::endl;
    } else {
        std::cout << "=====================" << std::endl;
    }


    return samples;
}

