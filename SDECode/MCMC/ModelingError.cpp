#include "ModelingError.hpp"

ModErr::ModErr(oneDimSde sdeCoarse, oneDimSde sdeFine, double (*V1) (double), double IC,
               VectorXd& priorMean, VectorXd& priorStdDev, double T,
               unsigned int N, std::vector<double>& observations, double noise, bool homogen):
        sdeCoarse(sdeCoarse),
        sdeFine(sdeFine),
        V1(V1),
        IC(IC),
        priorMean(priorMean),
        priorStdDev(priorStdDev),
        T(T),
        N(N),
        obs(observations),
        noise(noise),
        homogen(homogen)
{}

void ModErr::computePF(unsigned int nParam, unsigned int nMC)
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
    // #pragma omp parallel for num_threads(5)
    for (unsigned int i = 0; i < nParam; i++) {
        std::shared_ptr<ForwardPFModErr> FwdModErr;
        FwdModErr = std::make_shared<ForwardPFModErr>(sdeFine, sdeCoarse, h, seedMC, V1, homogen);
        ParFilMod PFMod(obs, T, IC, noise, param(0), nMC, FwdModErr);

        if (nParam == 1) {
            for (unsigned int j = 0; j < param.size(); j++) {
                param(j) = priorMean(j);
            }
        } else {
            for (unsigned int j = 0; j < param.size(); j++) {
                param(j) = priorMean(j) + priorStdDev(j) * gaussian(seedParam);
            }
        }

        PFMod.compute(param, true);
        auto solAndModeling = PFMod.sampleX();

        for (unsigned int k = 0; k < N + 1; k++) {
            errors[i][k] = solAndModeling[k](1);
        }
    }

    weights.resize(nParam);
    for (unsigned int i = 0; i < nParam; i++) {
        weights[i] = 1.0 / nParam;
    }

}

void ModErr::computePFAlt(unsigned int nMC)
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
    FwdModErr = std::make_shared<ForwardPFModErr>(sdeFine, sdeCoarse, h, seedMC, V1, homogen);
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

void ModErr::getModErr(std::vector<std::vector<double>>& getData)
{
    getData = errors;
}

void ModErr::getWeights(std::vector<double>& getData)
{
    getData = weights;
}