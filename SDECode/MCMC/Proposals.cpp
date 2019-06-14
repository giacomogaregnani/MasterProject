#include "Proposals.hpp"

Proposals::Proposals(double stdDev, std::vector<double> factors):
    proposalStdDev(stdDev),
    stdDev(stdDev),
    factors(factors)
{
    proposalNormal = std::normal_distribution<double>(0.0, 1.0);
}

VectorXd Proposals::RWProposal(VectorXd& theta, std::default_random_engine* generator)
{
    long int nParam = theta.size();
    VectorXd proposedValue(nParam);

    if (factors.empty()) {
        factors = std::vector<double>(nParam, 1.0);
    }

    for (long int i = 0; i < nParam; i++) {
        proposedValue(i) = theta(i) + factors[i] * stdDev * proposalNormal(*generator);
    }

    return proposedValue;
}

VectorXd Proposals::genSample(VectorXd &theta, std::default_random_engine *generator)
{
    return RWProposal(theta, generator);
}
