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

    double T = 100;
    auto N = static_cast<unsigned int>(pow(2, 16));

    std::vector<double> epsTestSet = {0.04};
    for (unsigned int i = 0; i < 8; i++)
        epsTestSet.push_back(epsTestSet[i] + 0.02);

    std::ofstream output(DATA_PATH + std::string("testParam.txt"), std::ofstream::out | std::ofstream::trunc);
    std::ofstream outputAvg(DATA_PATH + std::string("testParamAvg.txt"), std::ofstream::out | std::ofstream::trunc);
    std::ofstream outputSol(DATA_PATH + std::string("testParamSol.txt"), std::ofstream::out | std::ofstream::trunc);
    std::ofstream outputTest(DATA_PATH + std::string("testParamTest.txt"), std::ofstream::out | std::ofstream::trunc);

    // Compute coefficients of the homogenised equation
    std::vector<double> homCoeffs = computeHomCoeffs(param, (2.0*M_PI), &V1);
    output << homCoeffs[0] << "\t" << homCoeffs[1] << std::endl;
    outputAvg << homCoeffs[0] << "\t" << homCoeffs[1] << std::endl;

    // Compute the estimators for different value of epsilon
    double sigmaHat, AHat;

    unsigned int deltaRatio = 2;
    unsigned int nDeltas = 8;

    std::vector<double> avg;

    for (auto it : epsTestSet) {
        param(0) = it;

        std::cout << "============" << std::endl;
        std::cout << "\u03B5 = " << it << std::endl;
        output << param(0) << "\t";
        outputAvg << param(0) << "\t";

        auto solution = generateObservations1D(sde, 0.0, param, T, N, 0);
        for (const auto &itSol : solution)
            outputSol << std::fixed << std::setprecision(10) << itSol << "\t";
        outputSol << std::endl;

        auto NSampling = N;
        for (unsigned int j = 0; j < nDeltas; j++) {
            std::cout << "\u03B4 = " << T / NSampling << std::endl;

            // Estimate with sub-sampling
            std::vector<double> samples = {};
            auto increment = static_cast<unsigned int>(std::round(pow(deltaRatio, j)));
            for (unsigned int k = 0; k < N; k += increment)
                samples.push_back(solution[k]);

            sigmaHat = estimateSigma(samples, T / NSampling);
            AHat = estimateA(samples, T / NSampling, &gradV0, param);
            output << AHat << "\t" << sigmaHat << "\t";

            // Estimate with averaging
            avg = averageSequence(solution, increment);
            sigmaHat = estimateSigma(avg, T / N);
            AHat = estimateA(avg, T / N, &gradV0, param);
            outputAvg << AHat << "\t" << sigmaHat << "\t";

            NSampling = NSampling / deltaRatio;

            if (j == 2) {
                for (const auto &itSol : avg)
                    outputTest << std::fixed << std::setprecision(10) << itSol << "\t";
                outputTest << std::endl;
            }
        }
        output << std::endl;
        outputAvg << std::endl;
    }

    output.close();
    outputAvg.close();
    outputSol.close();
    outputTest.close();
    return 0;
}