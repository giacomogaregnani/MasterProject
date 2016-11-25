#include <cmath>
#include <vector>
#include <random>
#include <Solver.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include "problems.hpp"

MatrixXd diffusion(VectorXd x, std::vector<double>& param, double sigma, double h)
{
    long int size = x.size();
    MatrixXd M = MatrixXd::Identity(size, size) * h * sigma;
    return M;
}

MatrixXd func(MatrixXd a, MatrixXd b, MatrixXd I, double diffusion)
{
    return (a * b.transpose() + b * a + I * diffusion);
}

MatrixXd sqrtFunc(MatrixXd W, MatrixXd Fx, MatrixXd I, double diffusion)
{
    MatrixXd WInv = W.inverse();
    return Fx * W + diffusion * I * WInv.transpose();
}

int main(int argc, char* argv[]) {
    // =========================
    // INITIALIZATION
    // =========================
    problems problem = VDPOL;
    odeDef testODE;
    testODE.ode = problem;
    setProblem(&testODE);
    // =========================

    // =========================
    // Generate a "Gauss" trajectory
    // =========================
    std::string folder = DATA_PATH;
    std::string slash("/");
    std::ofstream meanResults;
    std::ofstream varResults;
    std::string meanPath = "meansHires.txt";
    std::string fullPath = folder + slash + meanPath;
    meanResults.open(fullPath, std::ofstream::out | std::ofstream::trunc);
    std::string varPath = "varHires.txt";
    fullPath = folder + slash + varPath;
    varResults.open(fullPath, std::ofstream::out | std::ofstream::trunc);

    MatrixXd initialVar = MatrixXd::Zero(testODE.size, testODE.size);
    // MatrixXd initialVar = 0.01 * MatrixXd::Identity(testODE.size, testODE.size);

    double h = 0.001, finalTime = 50.0, sigma = 1.0;
    int nSteps = static_cast<int> (finalTime / h);

    VectorXd mean = testODE.initialCond;
    MatrixXd variance = initialVar;
    MatrixXd W = initialVar;
    MatrixXd Fx(testODE.size, testODE.size);
    MatrixXd I = MatrixXd::Identity(testODE.size, testODE.size);
    double diffusion = h * h * sigma;

    meanResults << mean.transpose() << "\n";
    varResults << variance << "\n";

    VectorXd kOldM(testODE.size);
    VectorXd kNewM(testODE.size);
    VectorXd kCurrM(testODE.size);
    MatrixXd kOldV(testODE.size, testODE.size);
    MatrixXd kNewV(testODE.size, testODE.size);
    MatrixXd kCurrV(testODE.size, testODE.size);
    int nStages, minStages = 2;
    double coeff, coeff2;

    Fx = testODE.odeJac(mean, testODE.refParam);
    double rho;

    for (int i = 0; i < nSteps; i++) {

        rho = 2.5 * powerMethod(Fx, 0.1, 150);
        // varResults << Fx << "\n";
        nStages = static_cast<int>(std::ceil(sqrt(0.5 * h * rho))) + 1;
        nStages = std::max(minStages, nStages);
        std::cout << rho / 2.5 << " " << nStages << std::endl;

        coeff = h / (nStages * nStages);
        coeff2 = 2.0 * coeff;

        kOldV = variance;
        // dkOldV = W;
        kOldM = mean;
        kNewV = variance + func(kOldV, Fx, I, diffusion) * coeff;
        // kNewV = W + sqrtFunc(W, Fx, I, diffusion) * coeff;
        kNewM = mean + testODE.odeFunc(mean, testODE.refParam) * coeff;

        for (int j = 2; j < nStages + 1; j++) {
            kCurrM = testODE.odeFunc(kNewM, testODE.refParam) * coeff2 +
                     kNewM * 2.0 - kOldM;
            Fx = testODE.odeJac(kNewM, testODE.refParam);
            kCurrV = func(kNewV, Fx, I, diffusion) * coeff2 +
                     kNewV * 2.0 - kOldV;
            // kCurrV = sqrtFunc(W, Fx, I, diffusion) * coeff2 +
            //      kNewV * 2.0 - kOldV;
            kOldV = kNewV;
            kNewV = kCurrV;
            kOldM = kNewM;
            kNewM = kCurrM;
        }

        // W = kCurrV;
        //variance = W * W.transpose();
        mean = kCurrM;
        variance = kCurrV;

        meanResults << mean.transpose() << "\n";
        varResults << variance << "\n";
    }

    meanResults.close();
    varResults.close();
    return 0;
}