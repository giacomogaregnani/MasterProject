#include <cmath>
#include <random>
#include <iostream>
#include "mcmcTools.hpp"
#include <Solver.hpp>
#include <Eigen/Cholesky>

#define PI 3.1415926535897
#define M_EPS 0.0000001

std::vector<VectorXd> metropolisHastings(VectorXd initialCond, std::vector<double>& param, double sigma,
                                         int size, double h, double finalTime,
                                         std::vector<VectorXd>& data,
                                         std::vector<double>& dataTimes,
                                         VectorXd (*odeFunc) (VectorXd, std::vector<double>&),
                                         VectorXd priorMean, VectorXd priorVariance, int internalMC,
                                         int nStepsMC, double varData, double* accRatio)
{
    *accRatio = 0.0;
    std::default_random_engine wGenerator{(unsigned int) time(NULL)};
    std::default_random_engine solGenerator{(unsigned int) time(NULL)};
    std::normal_distribution<double> normal(0.0, 1.0);
    std::uniform_real_distribution<double> unif;
    size_t nData = data.size();

    double oldLike = 0.0;
    int nParam = (int) param.size();
    VectorXd paramVec(nParam);
    for (int i = 0; i < nParam; i++) {
        paramVec(i) = param[i];
    }
    double l;
    double prior;

    // Evaluate posterior for initial guess
    ProbMethod<EulerForward> firstSolver(size, h, initialCond, param, odeFunc, sigma);
    for (int j = 0; j < internalMC; j++) {
        double t = 0;
        int count = 0;
        std::vector<VectorXd> solutionAtDataTimes(nData, VectorXd(size));
        while (t < finalTime) {
            firstSolver.oneStep(solGenerator, h);
            t = t + h;
            if (std::abs(t - dataTimes[count]) < h / 10.0) {
                std::cout << t << std::endl;
                solutionAtDataTimes[count] = firstSolver.getSolution();
                count++;
            }
        }
        oldLike += evalLogLikelihood(data, solutionAtDataTimes, size, varData);
    }
    oldLike /= internalMC;
    double oldPrior = evalLogPrior(paramVec, priorMean, priorVariance, nParam);

    VectorXd paramVecN(nParam);
    std::vector<VectorXd> mcmcPath(nStepsMC, VectorXd(nParam));
    std::vector<VectorXd> solutionAtDataTimes(nData, VectorXd(size));
    VectorXd w(nParam);
    std::vector<double> paramStd(nParam);

    // =========== RAM UPDATE ===========
    double gamma = 0.01;
    MatrixXd I = MatrixXd::Identity(nParam, nParam);
    MatrixXd initialCov = gamma * I;
    LLT<MatrixXd> chol(initialCov);
    MatrixXd S = chol.matrixL();
    double desiredAlpha = 0.25;
    MatrixXd tmp(nParam, nParam);
    // ==================================

    for (int i = 0; i < nStepsMC; i++)
    {
        if (i % 50 == 0) {
            std::cout << "iteration: " << i << std::endl;
            std::cout << "likelihood: " << oldLike << std::endl;
            std::cout << "S matrix: " << std::endl << S << std::endl;
            std::cout << "prior: " << oldPrior << std::endl;
            std::cout << "param: " << paramVec.transpose() << std::endl;
        }

        for (int j = 0; j < nParam; j++) {
            w(j) = normal(wGenerator);
        }
        paramVecN = paramVec + S * w;

        for (int j = 0; j < nParam; j++) {
            paramStd[j] = paramVecN(j);
        }

        double lVec[internalMC];
        int MCindex;
#pragma omp parallel for num_threads(20) private(MCindex)
        for (MCindex = 0; MCindex < internalMC; MCindex++) {
            std::vector<VectorXd> solutionAtDataTimesPar(nData, VectorXd(size));
            ProbMethod<EulerForward> solver(size, h, initialCond, paramStd, odeFunc, sigma);
            double t = 0;
            int count = 0;
            while (t < finalTime) {
                solver.oneStep(solGenerator, h);
                t = t + h;
                if (std::abs(t - dataTimes[count]) < h / 10.0) {
                    solutionAtDataTimesPar[count] = solver.getSolution();
                    count++;
                }
            }
            lVec[MCindex] = evalLogLikelihood(data, solutionAtDataTimesPar, size, varData);
        }
        l = 0.0;
        for (int j = 0; j < internalMC; j++) {
            l += lVec[j];
        }
        l /= internalMC;

        // Add computations with old value of parameters
        for (int j = 0; j < nParam; j++) {
            paramStd[j] = paramVec(j);
        }

        double oldParamL = 0.0;
#pragma omp parallel for num_threads(20) private(MCindex)
        for (MCindex = 0; MCindex < internalMC; MCindex++) {
            std::vector<VectorXd> solutionAtDataTimesPar(nData, VectorXd(size));
            ProbMethod<EulerForward> solverOld(size, h, initialCond, paramStd, odeFunc, sigma);
            double t = 0;
            int count = 0;
            while (t < finalTime) {
                solverOld.oneStep(solGenerator, h);
                t = t + h;
                if (std::abs(t - dataTimes[count]) < h / 10.0) {
                    solutionAtDataTimesPar[count] = solverOld.getSolution();
                    count++;
                }
            }
            lVec[MCindex] = evalLogLikelihood(data, solutionAtDataTimesPar, size, varData);
        }

        for (int j = 0; j < internalMC; j++) {
            oldParamL += lVec[j];
        }
        oldParamL /= internalMC;

        prior = evalLogPrior(paramVecN, priorMean, priorVariance, nParam);
        double u = unif(wGenerator);
        double alpha = std::min<double>(1.0, exp((l + prior) - (oldParamL + oldPrior)));
        if (u < alpha) {
            *accRatio += 1.0 / nStepsMC;
            paramVec = paramVecN;
            oldLike = l;
            oldPrior = prior;
            gamma = std::min<double>(gamma * 2.0, 1.0);
        } else {
            gamma = std::max<double>(gamma / 2.0, 0.01);
        }
        mcmcPath[i] = paramVec;

        // RAM UPDATE
        MatrixXd WWT = w * w.transpose();
        double WTW = w.dot(w);
        double gammaI = std::min(1.0, 2.0 * pow(static_cast<double>(i + 1), -0.75));
        double diffAlpha = alpha - desiredAlpha;
        double coeff = gammaI * diffAlpha / WTW;
        MatrixXd C = I + coeff * WWT;
        tmp = S * C * S.transpose();
        LLT<MatrixXd> cholIt(tmp);
        S = cholIt.matrixL();
    }
    return mcmcPath;
}
