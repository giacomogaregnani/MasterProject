#include "mcmcTools.hpp"
#include <cmath>
#include <random>
#include <iostream>

void StoppedMetropolis(std::vector<VectorXd>& mcmcPath,
                       odeDef odeModel, std::vector<double>& param, double sigma,
                       double h, double finalTime, std::vector<VectorXd>& data, std::vector<double>& dataTimes,
                       VectorXd priorMean, VectorXd priorVariance, int internalMC,
                       int nStepsMC, double varData, bool posPar,
                       std::vector<double>& likelihoods, std::default_random_engine& generator,
                       MatrixXd& S, int index)
{
    // Initialize stuff
    size_t nParam = param.size();
    size_t nData = data.size();
    VectorXd paramVecN(nParam);
    likelihoods.resize(static_cast<size_t>(nStepsMC));
    VectorXd w(nParam);
    std::vector<double> paramStd(nParam);
    std::normal_distribution<double> normal(0.0, 1.0);
    std::uniform_real_distribution<double> unif;

    int MCindex;
    double lVec[internalMC];
    double l, prior, oldPrior, oldLike;

    double desiredAlpha = 0.25;

    VectorXd paramVec(nParam);
    for (int i = 0; i < nParam; i++) {
        paramVec(i) = param[i];
    }

    // Compute the posterior on the initial guess
    ProbMethod<RungeKutta> solver(odeModel, h, param, sigma);
    std::vector<VectorXd> solutionAtDataTimes(nData, VectorXd(odeModel.size));
    for (MCindex = 0; MCindex < internalMC; MCindex++) {
        double t = 0;
        int count = 0;
        while (t < finalTime) {
            solver.oneStep(generator, h);
            t = t + h;
            if (std::abs(t - dataTimes[count]) < h / 10.0) {
                solutionAtDataTimes[count++] = solver.getSolution();
            }
        }
        oldLike += evalLogLikelihood(data, solutionAtDataTimes, odeModel.size, varData);
        solver.resetIC();
    }
    oldLike /= internalMC;
    oldPrior = evalLogPrior(paramVec, priorMean, priorVariance, nParam);

    // Only for positive parameters
    double oldGauss, newGauss;

    printf("beginning index : %d\n", index);

    for (int i = 1; i < nStepsMC; i++) {
        if (i % 50 == 0) {
            printf("iteration: %d \n", i);
            /* std::cout << "iteration: " << i << std::endl;
            std::cout << "likelihood: " << oldLike << std::endl;
            std::cout << "step: " << std::endl << S << std::endl;
            std::cout << "prior: " << oldPrior << std::endl;
            std::cout << "param: " << paramVec.transpose() << std::endl; */
        }

        // Generate new guess for parameter
        bool isPositive = false;
        while (!isPositive) {
            for (int j = 0; j < nParam; j++) {
                w(j) = normal(generator);
            }
            paramVecN = paramVec + S * w;
            if (posPar) {
                if (paramVecN(0) > 0) {
                    isPositive = true;
                    oldGauss = phi(paramVec(0) / (S(0, 0) * S(0, 0)));
                    newGauss = phi(paramVecN(0) / (S(0, 0) * S(0, 0)));
                }
            } else {
                isPositive = true;
            }
        }
        for (int j = 0; j < nParam; j++) {
            paramStd[j] = paramVecN(j);
        }


        // Compute likelihood for the new parameter
        for (MCindex = 0; MCindex < internalMC; MCindex++) {
            ProbMethod<RungeKutta> solverNew(odeModel, h, paramStd, sigma);
            std::vector<VectorXd> solutionAtDataTimesPar(nData, VectorXd(odeModel.size));
            double t = 0;
            int count = 0;
            while (t < finalTime) {
                solverNew.oneStep(generator, h);
                t = t + h;
                if (std::abs(t - dataTimes[count]) < h / 10.0) {
                    solutionAtDataTimesPar[count] = solverNew.getSolution();
                    count++;
                }
            }
            lVec[MCindex] = evalLogLikelihood(data, solutionAtDataTimesPar, odeModel.size, varData);
            solverNew.resetIC();
        }
        l = 0.0;
        for (int j = 0; j < internalMC; j++) {
            l += lVec[j];
        }
        l /= internalMC;

        // Compute prior on the new parameter
        prior = evalLogPrior(paramVecN, priorMean, priorVariance, nParam);

        // Return on the old parameter
        for (int j = 0; j < nParam; j++) {
            paramStd[j] = paramVec(j);
        }

        // Compute likelihood for the old parameter
        for (MCindex = 0; MCindex < internalMC; MCindex++) {
            ProbMethod<RungeKutta> solverOld(odeModel, h, paramStd, sigma);
            std::vector<VectorXd> solutionAtDataTimesPar(nData, VectorXd(odeModel.size));
            double t = 0;
            int count = 0;
            while (t < finalTime) {
                solverOld.oneStep(generator, h);
                t = t + h;
                if (std::abs(t - dataTimes[count]) < h / 10.0) {
                    solutionAtDataTimesPar[count] = solverOld.getSolution();
                    count++;
                }
            }
            lVec[MCindex] = evalLogLikelihood(data, solutionAtDataTimesPar, odeModel.size, varData);
            solverOld.resetIC();
        }
        double oldParamL = 0.0;
        for (int j = 0; j < internalMC; j++) {
            oldParamL += lVec[j];
        }
        oldParamL /= internalMC;

        // Generate (log) of uniform random variable
        double u = log(unif(generator));

        // Compute probability of acceptance
        double alpha;
        if (posPar) {
            alpha = std::min<double>(0.0, (l + prior + log(oldGauss)) - (oldParamL + oldPrior + log(newGauss)));
        } else {
            alpha = std::min<double>(0.0, (l + prior) - (oldParamL + oldPrior));
        }
        // Update the chain
        if (u < alpha) {
            paramVec = paramVecN;
            oldLike = l;
            oldPrior = prior;
        }

        mcmcPath.push_back(paramVec);

        // Update RAM
        S = RAMupdate(S, w, exp(alpha), desiredAlpha, nParam, index++);
    }
}