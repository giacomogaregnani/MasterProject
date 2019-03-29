#ifndef PROPOSALS_HPP
#define PROPOSALS_HPP

#include <Eigen/Dense>
#include <random>

using namespace Eigen;

class Proposals {
private:

    std::normal_distribution<double> proposalNormal;

    double proposalStdDev;

    std::vector<double> proposalParam;

    bool isCN;

public:

    Proposals() {};

    Proposals(double stdDev, std::vector<double> param = {});

    VectorXd RWProposal(VectorXd& theta, std::default_random_engine* generator);

    VectorXd CNProposal(VectorXd& theta, std::default_random_engine* generator);

    VectorXd genSample(VectorXd& theta, std::default_random_engine* generator);
};

#endif