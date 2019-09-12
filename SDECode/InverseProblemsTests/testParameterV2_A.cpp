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
    double IC = 10.0;

    VectorXd param(3);
    param(0) = 0.005;
    param(1) = 1.0; // True multiscale alpha
    param(2) = 1.0; // True multiscale sigma
    double eps = param(0);

    std::ofstream output(DATA_PATH + std::string("testParamZeta.txt"), std::ofstream::out | std::ofstream::trunc);

    // Compute coefficients of the homogenised equation
    std::vector<double> homCoeffs = computeHomCoeffs(param, (2.0*M_PI), &V1);
    output << homCoeffs[0] << "\t" << homCoeffs[1] << std::endl;

    // Compute the estimators for different value of epsilon
    double AHat;
    std::vector<double> avg;

    // Refer to the Caltech notes
    double beta = 2.0;
    VectorXd zetaVec;
    zetaVec.setLinSpaced(21, beta-1.0, beta);

    for (unsigned int i = 0; i < zetaVec.size(); i++) {
        double zeta = zetaVec(i);
        double gamma = 2.0*beta - zeta;

        auto h = std::pow(eps, beta);
        auto N = static_cast<unsigned int>(std::round(std::pow(eps, -gamma)));
        auto T = h * N;
        auto delta = static_cast<unsigned int>(std::round(std::pow(eps, -zeta)));

        std::cout << "============" << std::endl;
        std::cout << zeta << std::endl;
        output << zeta << "\t";

        auto solution = generateObservations1D(sde, IC, param, T, N, 0);
        avg = averageSequence(solution, delta);
        AHat = estimateA(avg, h, &gradV0);
        output << AHat << std::endl;
        std::cout << "A = " << AHat << std::endl;
    }

    output.close();
    return 0;
}