#ifndef PROPOSALS_HPP
#define PROPOSALS_HPP

#include <Eigen/Dense>
#include <random>

using namespace Eigen;

class Proposals {
private:
    std::normal_distribution<double> proposalNormal;
    double proposalStdDev;
    bool isRAM;
    double alphaStar;
    double stdDev;
    MatrixXd stdDevRAM;

public:
    Proposals() = default;
    Proposals(double stdDev, double alphaStar = 0);
    VectorXd RWProposal(VectorXd& theta, std::default_random_engine* generator);
    VectorXd RAMProposal(VectorXd& theta, std::default_random_engine* generator);
    VectorXd genSample(VectorXd& theta, std::default_random_engine* generator);
};

#endif