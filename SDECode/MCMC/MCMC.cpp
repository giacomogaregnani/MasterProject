#include "MCMC.hpp"

MCMC::MCMC(VectorXd& initGuess,
           std::shared_ptr<Proposals> proposal,
           std::shared_ptr<Posterior> posterior,
           unsigned long nMCMC):
        proposal(proposal),
        posterior(posterior),
        initGuess(initGuess),
        nMCMC(nMCMC)
{
    sizeParam = initGuess.size();
    samples.resize(nMCMC+1);
    samples[0] = initGuess;
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
    unsigned long accRatio = 0;

    VectorXd proposedValue(sizeParam), currValue(sizeParam), currEst(sizeParam);

    // Compute posterior on first guess
    oldPosterior = posterior->computePosterior(samples[0]);

    currEst = initGuess;

    for (unsigned long i = 0; i < nMCMC; i++) {
        currValue = samples[i];

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
            samples[i+1] = proposedValue;
            oldPosterior = newPosterior;
            ++accRatio;
        } else {
            samples[i+1] = currValue;
        }
        currEst += samples[i+1];

        if (i % 100 == 0 && verbose) {
            std::cout << "Completed " << i << " iterations out of " << nMCMC << std::endl
                      << "Accepted (log)posterior = " << oldPosterior << std::endl
                      << "Current acc. ratio = " << 1.0 / (i + 1) * accRatio << std::endl
                      << "Current estimate = " << currEst.transpose() / (i + 2) << std::endl;
        }
    }

    // Take 10% burn-in out
    unsigned long burnIn = 0; // samples.size() / 10;
    samples.erase(samples.begin(), samples.begin() + burnIn);

    // Compute MCMC estimate
    VectorXd finalEstimate = VectorXd::Zero(sizeParam);
    for (auto const &it : samples) {
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