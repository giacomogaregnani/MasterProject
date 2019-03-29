#include "Proposals.hpp"

Proposals::Proposals(double stdDev, std::vector<double> param):
    proposalStdDev(stdDev),
    proposalParam(param)
{
    proposalNormal = std::normal_distribution<double>(0.0, stdDev);
    isCN = false;
    if (!param.empty()){
        isCN = true;
    }
}

VectorXd Proposals::RWProposal(VectorXd& theta, std::default_random_engine* generator)
{
    long int nParam = theta.size();
    VectorXd proposedValue(nParam);

    for (long int i = 0; i < nParam; i++) {
        proposedValue(i) = theta(i) + proposalNormal(*generator);
    }

    return proposedValue;
}

VectorXd Proposals::CNProposal(VectorXd &theta, std::default_random_engine *generator)
{
    long int nParam = theta.size();
    VectorXd proposedValue(nParam);

    for (long int i = 0; i < nParam; i++) {
        proposedValue(i) = sqrt(1 - proposalParam[0] * proposalParam[0]) * theta(i)
                           + proposalParam[0] * proposalNormal(*generator);
    }

    return proposedValue;
}

VectorXd Proposals::genSample(VectorXd &theta, std::default_random_engine *generator)
{
    if (isCN)
        return CNProposal(theta, generator);
    return RWProposal(theta, generator);
}
