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

double gradV1(double x)
{
    return -std::sin(x);
}

double multiDrift(double x, VectorXd& p)
{
    return -1.0 * (std::exp(p(1)) * gradV0(x) + 1 / p(0) * gradV1(x / p(0)));
}

double homoDrift(double x, VectorXd& p)
{
    return -1.0 * std::exp(p(1)) * gradV0(x);
}

double diffusion(double x, VectorXd &p)
{
    // return std::sqrt(2.0 * 0.5);
    return std::sqrt(2.0 * std::exp(p(2)));
}

int main(int argc, char* argv[])
{
    oneDimSde sde{&multiDrift, &diffusion};
    oneDimSde sdeHomo{&homoDrift, &diffusion};

    VectorXd param(3);
    param(0) = 0.1; // epsilon
    param(1) = std::log(1.0); // True multiscale alpha
    param(2) = std::log(0.5); // True multiscale betainv

    double T = 1;
    unsigned int N = 1000;

    std::ofstream output(DATA_PATH + std::string("test2.txt"), std::ofstream::out | std::ofstream::trunc);
    std::ofstream outputSol(DATA_PATH + std::string("testSol2.txt"), std::ofstream::out | std::ofstream::trunc);
    std::ofstream outputLik(DATA_PATH + std::string("testLik2.txt"), std::ofstream::out | std::ofstream::trunc);

    unsigned long M = 100;
    double noise = 1e-4;
    double IC = 1.0;
    std::default_random_engine noiseSeed{2018};
    std::normal_distribution<double> noiseDistribution(0.0, noise);

    // Generate and perturb observations
    auto x = generateObservations1D(sdeHomo, IC, param, T, N, 0);
    outputSol << x[0] << "\t";
    for (unsigned long i = 1; i < x.size(); i++) {
        x[i] += noiseDistribution(noiseSeed);
        outputSol << x[i] << "\t";
    }
    outputSol << std::endl;

    std::shared_ptr<ParFil> PF;
    PF = std::make_shared<ParFil>(x, T, IC, 1, noise, sdeHomo, param(0), M);
    std::vector<std::vector<double>> X;

    for (unsigned int i = 0; i < 1; i++) {
        PF->compute(param);
        outputLik << PF->getLikelihood() << std::endl;
        if ((i+1) % 100 == 0) {
            std::cout << "iteration " << i+1 << std::endl;
        }
    }

    X = PF->getX();
    for (auto const &it : X) {
        for (auto const &itit : it)
            output << itit << "\t";
        output << std::endl;
    }
    for (auto const &it : PF->getBestX()) {
        output << it << "\t";
    }
    output << std::endl;

    output.close();
    outputSol.close();
    outputLik.close();
    return 0;
}