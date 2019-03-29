#include "MCMC.hpp"

MCMC::MCMC(VectorXd& initGuess,
           Proposals* proposal,
           Posterior* posterior,
           unsigned long nMCMC):
        proposal(proposal),
        posterior(posterior),
        nMCMC(nMCMC),
        initGuess(initGuess)
{
    samples = {initGuess};
    probGen = std::uniform_real_distribution<double>(0.0, 1.0);
}

void MCMC::eraseSample(void)
{
    samples.clear();
    samples.push_back(initGuess);
}

std::vector<VectorXd>& MCMC::compute(std::default_random_engine* proposalSeed,
                                     std::default_random_engine* acceptanceSeed,
                                     bool noisy, bool verbose)
{
    double newPosterior, oldPosterior, alpha;
    VectorXd currEstimate;
    currEstimate = samples.back();
    unsigned long accRatio = 0;

    VectorXd proposedValue(sizeParam), currValue(sizeParam);

    // Compute posterior on first guess
    oldPosterior = posterior->computePosterior(samples.back());

    for (unsigned long i = 0; i < nMCMC; i++) {

        currValue = samples.back();

        // Generate new guess
        proposedValue = proposal->genSample(currValue, proposalSeed);

        // Evaluate posterior on new guess
        newPosterior = posterior->computePosterior(proposedValue);

        // For the noisy algorithm, re-evaluate on old guess
        if (noisy)
            oldPosterior = posterior->computePosterior(currValue);

        // compute probability
        alpha = newPosterior - oldPosterior;
        alpha = std::min(0.0, alpha);

        // Update chain
        if (std::log(probGen(*acceptanceSeed)) < alpha) {
            samples.push_back(proposedValue);
            oldPosterior = newPosterior;
            ++accRatio;
        } else {
            samples.push_back(currValue);
        }

        // Update current estimate
        currEstimate = currEstimate * static_cast<double>(i + 1) / (i + 2)
                       + samples.back() * 1.0 / (i + 2);

        if (i % 1000 == 0 && verbose) {
            std::cout << "Completed " << i << " iterations out of " << nMCMC << std::endl
                      // << "Current estimate = " << currEstimate.transpose() << std::endl
                      << "Accepted (log)posterior = " << oldPosterior << std::endl
                      << "Current acc. ratio = " << 1.0 / (i + 1) * accRatio << std::endl;
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
              << " ACCEPTANCE RATIO: " << static_cast<double>(accRatio) / nMCMC << std::endl
              << "=====================" << std::endl;

    return samples;
}