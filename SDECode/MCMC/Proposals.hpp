#ifndef PROPOSALS_HPP
#define PROPOSALS_HPP

#include <Eigen/Dense>
#include <random>

using namespace Eigen;

class Proposals {
private:
    std::normal_distribution<double> proposalNormal;
    double proposalStdDev;
    double stdDev;
    std::vector<double> factors;

public:
    Proposals() = default;
    Proposals(double stdDev, std::vector<double> factors = {});
    VectorXd RWProposal(VectorXd& theta, std::default_random_engine* generator);
    VectorXd genSample(VectorXd& theta, std::default_random_engine* generator);
};

#endif