#include <fstream>
#include <iostream>
#include <iomanip>
#include "generateObservations.hpp"
#include "computeHomogeneous.hpp"

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
    return std::sqrt(2.0 * std::exp(p(2)));
}

int main(int argc, char* argv[])
{
    oneDimSde sde;
    sde.drift = &multiDrift;
    sde.diffusion = &diffusion;
    oneDimSde sdeHomo;
    sdeHomo.drift = &homoDrift;
    sdeHomo.diffusion = &diffusion;

    double T = 10;
    unsigned int N = 8000;

    VectorXd tmpParam(3);
    tmpParam(0) = 0.04;  // Epsilon
    tmpParam(1) = 1.0;  // True multiscale alpha
    tmpParam(2) = 2.0;  // True multiscale betainv

    std::ofstream outputSol(DATA_PATH + std::string("TestHomog.txt"), std::ofstream::out | std::ofstream::trunc);

    // Compute coefficients of the homogenised equation
    std::vector<double> homCoeffs = computeHomCoeffs(tmpParam, (2.0*M_PI), &V1);

    // ============= Transform for positiveness ============= //
    VectorXd param = tmpParam;
    param(1) = std::log(param(1));
    param(2) = std::log(param(2));
    // ====================================================== //

    double IC = 0.0;
    std::random_device dev;
    unsigned int obsSeed;

    for (unsigned int M = 0; M < 10000; M++) {
        obsSeed = dev();
        auto x = generateObservations1D(sde, IC, param, T, N, obsSeed);
        outputSol << x[N] << "\t";

        // Compute a homogeneous path with the same Brownian for comparison
        VectorXd homoParam(3);
        homoParam(0) = 0.1;
        homoParam(1) = std::log(homCoeffs[0]);
        homoParam(2) = std::log(homCoeffs[1]);
        auto xHom = generateObservations1D(sdeHomo, IC, homoParam, T, N, obsSeed);
        outputSol << xHom[N] << "\t";
        outputSol << std::endl;
    }

    outputSol.close();
    return 0;
}