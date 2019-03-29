#include <fstream>
#include <iostream>
#include <ParFil.hpp>
#include "generateObservations.hpp"

// Remark: epsilon = p(0)

double gradV0(double x)
{
    return x;
}

double V1(double x)
{
    return std::cos(x);
}

double drift(double x, VectorXd& p)
{
    return -1.0 * p(1) * gradV0(x);
}

double diffusion(double x, VectorXd &p)
{
    return std::sqrt(2.0 * 0.5);
}

int main(int argc, char* argv[])
{
    oneDimSde sde{&drift, &diffusion};

    VectorXd param(3);
    param(0) = 0.0; // Fake
    param(1) = 1.0; // True multiscale alpha
    param(2) = 0.5; // True multiscale betainv

    double T = 1;
    unsigned int N = 1000;

    std::ofstream output(DATA_PATH + std::string("test.txt"), std::ofstream::out | std::ofstream::trunc);
    std::ofstream outputSol(DATA_PATH + std::string("testSol.txt"), std::ofstream::out | std::ofstream::trunc);

    unsigned long M = 500;
    double noise = 1e-2;
    double IC = 1.0;
    std::default_random_engine noiseSeed{(unsigned int) time(nullptr)};
    std::normal_distribution<double> noiseDistribution(0.0, noise);

    // Generate and perturb observations
    auto x = generateObservations1D(sde, IC, param, T, N);
    outputSol << x[0] << "\t";
    for (unsigned long i = 1; i < x.size(); i++) {
        x[i] += noiseDistribution(noiseSeed);
        outputSol << x[i] << "\t";
    }
    outputSol << std::endl;

    std::shared_ptr<ParFil> PF;
    PF = std::make_shared<ParFil>(x, T, IC, 1, noise, sde, 0.0, M);
    PF->compute(param);

    std::vector<double> bestGuess = PF->getBestX();

    for (auto const &it : bestGuess)
        output << it << std::endl;

    std::cout << "likelihood = " << PF->getLikelihood() << std::endl;
    output.close();
    return 0;
}