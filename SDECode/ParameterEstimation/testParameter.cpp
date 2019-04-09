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

    VectorXd param(3);
    param(1) = 1.0; // True multiscale alpha
    param(2) = 0.5; // True multiscale sigma

    double T = 1000;
    auto N = static_cast<unsigned int>(std::round(pow(2, 14)));

    std::vector<double> epsTestSet = {0.04};
    for (unsigned int i = 0; i < 8; i++) {
        epsTestSet.push_back(epsTestSet[i] + 0.02);
    }

    std::ofstream output(DATA_PATH + std::string("testParam.txt"), std::ofstream::out | std::ofstream::trunc);
    std::ofstream outputSol(DATA_PATH + std::string("testParamSol.txt"), std::ofstream::out | std::ofstream::trunc);

    // Compute coefficients of the homogenised equation
    std::vector<double> homCoeffs = computeHomCoeffs(param, (2.0*M_PIf64), &V1);
    output << 0.0 << "\t" << homCoeffs[0] << "\t" << homCoeffs[1] << std::endl;

    // Compute the estimators for different value of epsilon
    double sigmaHat, AHat, ABayes;

    for (auto it : epsTestSet) {
        param(0) = it;

        std::cout << "============" << std::endl;
        std::cout << "\u03B5 = " << it << std::endl;
        output << param(0) << "\t";

        auto solution = generateObservations1D(sde, 0.0, param, T, N);
        for (const auto &itSol : solution)
            outputSol << std::fixed << std::setprecision(10) << itSol << "\t";
        outputSol << std::endl;

        auto NSampling = N;

        for (unsigned int j = 0; j < 4; j++) {
            std::cout << "\u03B4 = " << T / NSampling << std::endl;
            std::vector<double> samples = {};
            for (unsigned int k = 0; k < N; k += static_cast<unsigned int>(std::round(pow(4, j))))
                samples.push_back(solution[k]);

            sigmaHat = estimateSigma(samples, T / NSampling);
            AHat = estimateA(samples, T / NSampling, &gradV0, param);
            ABayes = estimateABayes(samples, T / NSampling, &gradV0, param, param(1), 0.01, 1.0 / param(2));

            output << AHat << "\t" << sigmaHat << "\t" << ABayes << "\t";
            NSampling = NSampling / 4;
        }
        output << std::endl;
    }

    output.close();
    return 0;
}