#include "Proposals.hpp"

Proposals::Proposals(double stdDev, double alphaStar):
    proposalStdDev(stdDev),
    stdDev(stdDev),
    alphaStar(alphaStar)
{
    proposalNormal = std::normal_distribution<double>(0.0, 1.0);
    isRAM = false;
    if (alphaStar > 0){
        isRAM = true;
    }
}

VectorXd Proposals::RWProposal(VectorXd& theta, std::default_random_engine* generator)
{
    long int nParam = theta.size();
    VectorXd proposedValue(nParam);

    for (long int i = 0; i < nParam; i++) {
        proposedValue(i) = theta(i) + stdDev * proposalNormal(*generator);
    }

    return proposedValue;
}

VectorXd Proposals::RAMProposal(VectorXd& theta, std::default_random_engine* generator)
{
    long int nParam = theta.size();
    VectorXd proposedValue(nParam);


    return proposedValue;
}

VectorXd Proposals::genSample(VectorXd &theta, std::default_random_engine *generator)
{
    if (isRAM)
        return RAMProposal(theta, generator);
    return RWProposal(theta, generator);
}
