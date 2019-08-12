#include "ModelingErrorAll.hpp"
#include "../matplotlib-cpp-master/matplotlibcpp.h"

namespace plt = matplotlibcpp;

ModErrAll::ModErrAll(oneDimSde sdeCoarse, oneDimSde sdeFine, double (*V1) (double), double IC,
               VectorXd& priorMean, VectorXd& priorStdDev, double T,
               unsigned int N, std::vector<double>& observations, double noise):
        sdeCoarse(sdeCoarse),
        sdeFine(sdeFine),
        V1(V1),
        IC(IC),
        priorMean(priorMean),
        priorStdDev(priorStdDev),
        T(T),
        N(N),
        obs(observations),
        noise(noise)
{}

/*
void ModErrAll::computePF(unsigned int nParam, unsigned int nMC)
{
    errors.resize(nMC*nParam+1);
    for (unsigned int i = 0; i < errors.size(); i++) {
        errors[i].resize(N + 1);
        for (unsigned int j = 0; j < N + 1; j++) {
            errors[i][j] = 0.0;
        }
    }

    double BM;
    auto h = T / N;

    std::random_device dev;
    std::default_random_engine seedMC{dev()};
    std::default_random_engine seedParam{dev()};
    VectorXd param(priorMean.size());
    VectorXd paramHom(priorMean.size());
    VectorXd tmp;
    param(0) = priorMean(0);
    paramHom(0) = priorMean(0);

    std::normal_distribution<double> gaussian(0.0, 1.0);

    EM1D solverCoarse(sdeCoarse, param, seedMC);
    ParFil parFil(obs, T, IC, 1, noise, sdeFine, param(0), nMC);

    std::vector<double> timeVec(N + 1);
    for (unsigned int i = 0; i < N + 1; i++) {
        timeVec[i] = h * i;
    }

    // Compute the modeling error
    for (unsigned int i = 0; i < nParam; i++) {
        for (unsigned int j = 0; j < param.size(); j++) {
            param(j) = priorMean(j) + priorStdDev(j) * gaussian(seedParam);
        }
        tmp = computeHomogeneous(param, (2.0 * M_PI), V1);
        paramHom(1) = tmp(0);
        paramHom(2) = tmp(1);
        solverCoarse.modifyParam(paramHom);
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

        std::vector<double> solCoarseOld(nMC);
        unsigned int ancestor;
        for (unsigned int k = 0; k < N; k++) {
            for (unsigned int j = 0; j < nMC; j++) {
                solCoarseOld[j] = solCoarse[j][k];
            }
            for (unsigned int j = 0; j < nMC; j++) {
                BM = BMVec[j][k];
                ancestor = treePF[j][k];
                solCoarse[j][k] = solCoarseOld[ancestor];
                solCoarse[j][k + 1] = solverCoarse.oneStepGivenNoise(h, solCoarse[j][k], BM);
                // errors[i][k + 1] += (solFine[j][k+1] - solCoarse[j][k+1]) / nMC;
                errors[i*nMC+j][k + 1] += (solFine[j][k+1] - solCoarse[j][k+1]);
            }
        }
    }

    for (unsigned int i = 0; i < errors.size()-1; i++) {
        for (unsigned int j = 0; j < N+1; j++) {
            errors.back()[j] += errors[i][j] / (errors.size()-1);
        }
    }
} */

void ModErrAll::computePF(unsigned int nParam, unsigned int nMC)
{
    errors.resize(nParam);
    for (unsigned int i = 0; i < errors.size(); i++) {
        errors[i].resize(N + 1);
        for (unsigned int j = 0; j < N + 1; j++) {
            errors[i][j] = 0.0;
        }
    }

    auto h = T / N;

    std::random_device dev;
    std::default_random_engine seedMC{dev()};
    std::default_random_engine seedParam{dev()};
    VectorXd param(priorMean.size());
    VectorXd paramHom(priorMean.size());
    VectorXd tmp;
    param(0) = priorMean(0);
    paramHom(0) = priorMean(0);

    std::normal_distribution<double> gaussian(0.0, 1.0);

    // Compute the modeling error
    #pragma omp parallel for num_threads(5)
    for (unsigned int i = 0; i < nParam; i++) {
        std::shared_ptr<ForwardPFModErr> FwdModErr;
        FwdModErr = std::make_shared<ForwardPFModErr>(sdeFine, sdeCoarse, h, seedMC, V1);
        ParFilMod PFMod(obs, T, IC, noise, param(0), nMC, FwdModErr);

        for (unsigned int j = 0; j < param.size(); j++) {
            param(j) = priorMean(j) + priorStdDev(j) * gaussian(seedParam);
        }
        PFMod.compute(param, true);
        auto solAndModeling = PFMod.sampleX();

        for (unsigned int k = 0; k < N + 1; k++) {
            errors[i][k] = solAndModeling[k](1);
        }
    }
}

void ModErrAll::computePFAlt(unsigned int nMC)
{
    errors.resize(nMC);
    for (unsigned int i = 0; i < errors.size(); i++) {
        errors[i].resize(N + 1);
        for (unsigned int j = 0; j < N + 1; j++) {
            errors[i][j] = 0.0;
        }
    }

    auto h = T / N;

    std::random_device dev;
    std::default_random_engine seedMC{dev()};
    std::default_random_engine seedParam{dev()};
    VectorXd param(priorMean.size());
    VectorXd paramHom(priorMean.size());
    VectorXd tmp;
    param = priorMean;

    std::normal_distribution<double> gaussian(0.0, 1.0);

    // Compute the modeling error
    std::shared_ptr<ForwardPFModErr> FwdModErr;
    FwdModErr = std::make_shared<ForwardPFModErr>(sdeFine, sdeCoarse, h, seedMC, V1);
    ParFilMod PFMod(obs, T, IC, noise, param(0), nMC, FwdModErr);

    PFMod.compute(param, true);
    auto solAndModeling = PFMod.getX();

    for (unsigned int i = 0; i < nMC; i++) {
        for (unsigned int k = 0; k < N + 1; k++) {
            errors[i][k] = solAndModeling[i][k](1);
        }
    }
    weights = PFMod.getW();
}

void ModErrAll::getModErr(std::vector<std::vector<double>>& getData)
{
    getData = errors;
}

void ModErrAll::getWeights(std::vector<double>& getData)
{
    getData = weights;
}