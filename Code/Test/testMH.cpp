#include "mcmcTools.hpp"
#include <cmath>
#include <random>

double posterior(VectorXd x)
{
    double tmp = x(0) * x(0) - x(1);
    return exp(-10 * tmp * tmp - pow(x(1) - 0.25, 4));
}

std::vector<VectorXd> testMetropolis(VectorXd oldGuess, int nIter, double* accRatio,
                                     double gamma, bool RAM, double desiredAlpha)
{
    // Normal generators
    std::normal_distribution<double> normal(0.0, 1.0);
    std::default_random_engine generator{(unsigned int) time(NULL)};
    std::uniform_real_distribution<double> unif;

    // Initial posterior
    double oldPost = posterior(oldGuess);

    // Number of parameters (hard-coded to 2)
    int nParam = 2;

    // New guess
    VectorXd w(nParam);
    VectorXd newGuess(nParam);
    double newPost;

    // MCMC Path and initialization
    std::vector<VectorXd> MCMC;
    MCMC.push_back(oldGuess);

    // Probability
    double alpha, u;

    // Only for RAM
    MatrixXd S;
    if (RAM) {
        S = RAMinit(gamma, desiredAlpha, nParam);
    }

    for (int i = 0; i < nIter; i++) {
        // generate new guess
        for (int j = 0; j < nParam; j++) {
            w(j) = normal(generator);
        }
        if (RAM) {
            newGuess = oldGuess + S * w;
        }
        else {
            newGuess = oldGuess + gamma * w;
        }

        // evaluate posterior
        newPost = posterior(newGuess);

        // compute probability
        alpha = std::min(1.0, newPost / oldPost);
        u = unif(generator);

        // Update chain
        if (u < alpha) {
            oldGuess = newGuess;
            oldPost = newPost;
        }
        MCMC.push_back(oldGuess);

        // Only for RAM
        if (RAM) {
            S = RAMupdate(S, w, alpha, desiredAlpha, nParam, i + 1);
        }

    }

    return MCMC;
}