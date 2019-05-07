#include "ModelingError.hpp"
#include "../matplotlib-cpp-master/matplotlibcpp.h"

namespace plt = matplotlibcpp;

ModErr::ModErr(oneDimSde sdeCoarse, oneDimSde sdeFine, double IC,
               VectorXd& priorMean, VectorXd& priorStdDev, double T,
               unsigned int N, std::vector<double>& observations, double noise):
           sdeCoarse(sdeCoarse),
           sdeFine(sdeFine),
           IC(IC),
           priorMean(priorMean),
           priorStdDev(priorStdDev),
           T(T),
           N(N),
           obs(observations),
           noise(noise)
{
    means.resize(N+1);
    stdDevs.resize(N+1);
    for (unsigned int i = 0; i < N+1; i++) {
        means[i] = 0.0;
        stdDevs[i] = 0.0;
    }
}

void ModErr::compute(unsigned int nParam, unsigned int nMC)
{
    errors.resize(nParam);
    for (unsigned int i = 0; i < nParam; i++) {
        errors[i].resize(N);
        for (unsigned int j = 0; j < N; j++) {
            errors[i][j] = 0.0;
        }
    }

    double solFine, solCoarse, BM;
    auto h = T / N;
    auto sqrth = std::sqrt(h);

    std::random_device dev;
    std::default_random_engine seedParam{dev()};
    VectorXd param(priorMean.size());

    std::normal_distribution<double> gaussian(0.0, 1.0);

    EM1D solverFine(sdeFine, param, seedParam);
    EM1D solverCoarse(sdeCoarse, param, seedParam);

    // Compute the modeling error
    for (unsigned int i = 0; i < nParam; i++) {
        for (unsigned int j = 0; j < param.size(); j++) {
            param(j) = priorMean(j) + priorStdDev(j) * gaussian(seedParam);
        }
        solverFine.modifyParam(param);
        solverCoarse.modifyParam(param);

        for (unsigned int j = 0; j < nMC; j++) {
            std::default_random_engine seedMC{2019};
            gaussian.reset();
            solFine = IC;
            solCoarse = IC;
            for (unsigned int k = 0; k < N; k++) {
                BM = sqrth * gaussian(seedMC);
                solFine = solverFine.oneStepGivenNoise(h, solFine, BM);
                solCoarse = solverCoarse.oneStepGivenNoise(h, solCoarse, BM);
                errors[i][k] += (solFine - solCoarse) / nMC;
            }
        }
    }

    // Compute modeling error statistics
    for (unsigned int i = 0; i < nParam; i++) {
        for (unsigned int j = 0; j < N; j++) {
            means[j] += errors[i][j] / nParam;
        }
    }

    double temp;
    for (unsigned int i = 0; i < nParam; i++) {
        for (unsigned int j = 0; j < N; j++) {
            temp = errors[i][j] - means[j];
            stdDevs[j] += temp * temp / (nParam - 1);
        }
    }

    for (unsigned int j = 0; j < N; j++) {
        stdDevs[j] = std::sqrt(stdDevs[j]);
    }

    means.insert(means.begin(), 0.0);
    stdDevs.insert(stdDevs.begin(), 0.0);
}

void ModErr::computePF(unsigned int nParam, unsigned int nMC)
{
    errors.resize(nParam);
    for (unsigned int i = 0; i < nParam; i++) {
        errors[i].resize(N+1);
        for (unsigned int j = 0; j < N+1; j++) {
            errors[i][j] = 0.0;
        }
    }

    double BM;
    auto h = T / N;

    std::random_device dev;
    std::default_random_engine seedMC{dev()};
    std::default_random_engine seedParam{dev()};
    VectorXd param(priorMean.size());
    param(0) = priorMean(0);

    std::normal_distribution<double> gaussian(0.0, 1.0);

    EM1D solverCoarse(sdeCoarse, param, seedMC);
    ParFil parFil(obs, T, IC, 1, noise, sdeFine, param(0), nMC);

    std::vector<double> timeVec(N+1);
    for (unsigned int i = 0; i < N+1; i++) {
        timeVec[i] = h*i;
    }

    // Compute the modeling error
    for (unsigned int i = 0; i < nParam; i++) {
        for (unsigned int j = 0; j < param.size(); j++) {
            param(j) = priorMean(j) + priorStdDev(j) * gaussian(seedParam);
        }
        solverCoarse.modifyParam(param);

        parFil.compute(param);
        auto W = parFil.getW();
        auto BMVec = parFil.getBM();
        auto solFine = parFil.getX();
        auto treePF = parFil.getTree();
        std::vector<std::vector<double>> solCoarse(nMC);

        for (unsigned int j = 0; j < nMC; j++) {
            solCoarse[j].resize(N + 1);
            solCoarse[j][0] = IC;
        }

        unsigned int ancestor;
        for (unsigned int k = 0; k < N; k++) {
            for (unsigned int j = 0; j < nMC; j++) {
                BM = BMVec[j][k];
                ancestor = treePF[j][k];
                solCoarse[j][k] = solCoarse[ancestor][k];
                solCoarse[j][k+1] = solverCoarse.oneStepGivenNoise(h, solCoarse[j][k], BM);
                errors[i][k+1] += (solFine[j][k+1] - solCoarse[j][k+1]) / nMC;
            }
        }
        /* for (unsigned int j = 0; j < nMC; j++) {
            std::vector<double> BMTest(N + 1);
            BMTest[0] = 0.0;
            for (unsigned int ii = 0; ii < N; ii++) {
                BMTest[ii + 1] = BMTest[ii] + BMVec[j][ii];
            }
            plt::subplot(3, 1, 1);
            plt::plot(timeVec, BMTest);
            plt::subplot(3, 1, 2);
            plt::plot(timeVec, solCoarse[j]);
            plt::subplot(3, 1, 3);
            plt::plot(timeVec, solFine[j], "grey");
        }
        plt::subplot(3, 1, 3);
        plt::plot(timeVec, obs, "k");
        plt::show();

        plt::plot(W);
        plt::show(); */
    }

    // Compute modeling error statistics
    for (unsigned int i = 0; i < nParam; i++) {
        for (unsigned int j = 0; j < N+1; j++) {
            means[j] += errors[i][j] / nParam;
        }
    }
    double temp;
    for (unsigned int i = 0; i < nParam; i++) {
        for (unsigned int j = 0; j < N+1; j++) {
            temp = errors[i][j] - means[j];
            stdDevs[j] += temp * temp / (nParam - 1);
        }
    }

    for (unsigned int j = 0; j < N+1; j++) {
        stdDevs[j] = std::sqrt(stdDevs[j]);
    }
}

void ModErr::getStats(std::vector<double>& getMean, std::vector<double>& getStdDev)
{
    getMean = means;
    getStdDev = stdDevs;
}