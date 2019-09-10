#include <fstream>
#include <iostream>
#include <iomanip>
#include "generateObservations.hpp"
#include "computeHomogeneous.hpp"
#include "computeEstimators.hpp"

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
    return -1.0 * (p(1) * gradV0(x) + 1 / p(0) * gradV1(x / p(0)));
}

double diffusion(double x, VectorXd &p)
{
    return std::sqrt(2.0 * p(2));
}

int main(int argc, char* argv[])
{
    oneDimSde sde{&multiDrift, &diffusion};
    double IC = 0.0;

    VectorXd param(3);
    param(1) = 1.0; // True multiscale alpha
    param(2) = 0.5; // True multiscale sigma

    std::vector<double> epsTestSet = {0.01};
    for (unsigned int i = 0; i < 19; i++) {
        // epsTestSet.push_back(epsTestSet[i] * 1.2);
        epsTestSet.push_back(epsTestSet[i] + 0.01);
    }

    std::ofstream output(DATA_PATH + std::string("testParamAvg.txt"), std::ofstream::out | std::ofstream::trunc);
    std::ofstream outputSol(DATA_PATH + std::string("testParamSol.txt"), std::ofstream::out | std::ofstream::trunc);
    bool flag = true;

    // Compute coefficients of the homogenised equation
    std::vector<double> homCoeffs = computeHomCoeffs(param, (2.0*M_PI), &V1);
    output << homCoeffs[0] << "\t" << homCoeffs[1] << std::endl;

    // Compute the estimators for different value of epsilon
    double sigmaHat, AHat;
    std::vector<double> avg;

    // Refer to the Caltech notes
    double beta = 3.0;
    double zeta = 2.8;
    double gamma = 4.0;

    for (auto it : epsTestSet) {

        param(0) = it;

        auto h = std::pow(it, beta);
        auto N = static_cast<unsigned int>(std::round(std::pow(it, -gamma)));
        auto T = h * N;
        auto delta = static_cast<unsigned int>(std::round(std::pow(it, -zeta)));

        std::cout << "============" << std::endl;
        std::cout << "\u03B5 = " << it << std::endl;
        std::cout << "N = " << N << ", T = " << T << std::endl;
        std::cout << "\u03B4 = " << delta << std::endl;

        auto solution = generateObservations1D(sde, IC, param, T, N, 0);
        avg = averageSequence(solution, delta);
        sigmaHat = delta * estimateSigma(avg, h);
        AHat = estimateA(avg, h, &gradV0);

        output << param(0) << "\t" << AHat << "\t" << sigmaHat << std::endl;

        if (N < 100000 && flag) {
            for (auto const &itSol : solution)
                outputSol << itSol << "\t";
            outputSol << std::endl;
            for (auto const &itSol : avg)
                outputSol << itSol << "\t";
            outputSol << std::endl;
            flag = false;
        }
    }

    output.close();
    outputSol.close();
    return 0;
}