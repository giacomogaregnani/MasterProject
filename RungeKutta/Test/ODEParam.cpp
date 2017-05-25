#include <RandomTimeStep.hpp>
#include <MCMC.hpp>
#include <UtilitiesGG.hpp>

#ifndef PI
#define PI 3.14159265359
#endif

static std::vector<double> obsTimes;
struct IntegrationParam {
    odeDef ODE;
    Butcher tableau;
    double h;
    double p;
};
static IntegrationParam integrationParam;
static std::default_random_engine forwardGenerator{(unsigned int) time(NULL)};
static VectorXd priorParam;

double likelihood(VectorXd& x, std::vector<VectorXd>& obs)
{
    RungeKuttaRandomH method(&forwardGenerator, integrationParam.ODE,
                             EigVecToStdVec(x), integrationParam.tableau,
                             integrationParam.h, integrationParam.p + 0.5);

    VectorXd solution = integrationParam.ODE.initialCond;

    std::vector<VectorXd> forwardMapResults;
    unsigned long N = static_cast<unsigned long> (obsTimes.back() / integrationParam.h);

    unsigned int count = 0;
    for (unsigned long i = 0; i < N; i++) {
        solution = method.oneStep(solution);
        if (integrationParam.h * (i + 1) == obsTimes[count]) {
            forwardMapResults.push_back(solution);
            ++count;
        }
    }

    double l = 0.0;
    for (unsigned int j = 0; j < obs.size(); j++) {
        l -= 1.0 / (2.0 * 1e-2) * (forwardMapResults[j] - obs[j]).dot(forwardMapResults[j] - obs[j]);
    }

    return l;
}

double gaussPrior(VectorXd& x)
{
    return -0.5 * (x - priorParam).dot(x - priorParam);
}

double logPrior(VectorXd& x)
{
    bool negative = true;
    for (unsigned int i = 0; i < x.size(); i++) {
        negative = negative * (x(i) > 0);
    }
    if (!negative) {
        return -1e30;
    }

    VectorXd logarX(x.size());
    for (unsigned int i = 0; i < x.size(); i++) {
        logarX(i) = std::log(x(i));
    }

    VectorXd fact = logarX - priorParam;
    return - 1.0 * fact.dot(fact) / 2.0;
}


int main(int argc, char* argv[])
{
    // ODE
    integrationParam.ODE.ode = LORENZ;
    setProblem(&integrationParam.ODE);

    // Integration parameters
    integrationParam.h = 0.02;
    integrationParam.p = 1.5;

    integrationParam.tableau = Butcher(EULERFORWARD, EXPLICIT, 0);

    // Load observation(s)
    std::vector<VectorXd> observations;
    std::string refFilename = std::string(DATA_PATH) + argv[1] + ".txt";
    loadObservations(obsTimes, observations,
                     refFilename, integrationParam.ODE.size);

    // Perturb observation(s)
    for (unsigned int i = 0; i < obsTimes.size(); i++) {
        observations[i] += normalZeroMeanRandVec(integrationParam.ODE.size,
                                                 &forwardGenerator,
                                                 1e-1);
    }

    // Set prior parameters
    priorParam = StdVecToEig(integrationParam.ODE.refParam);

    // Proposal parameter
    std::normal_distribution<double>::param_type param(0.0, 0.1);

    // Init prior and likelihood
    ProbDensFunc priorObj;
    priorObj.PDF = &gaussPrior;
    LikelihoodFunc likObj;
    likObj.LIK = &likelihood;

    // MCMC Param
    unsigned long nMCMC = 50000;
    bool RAM = true;
    double desiredAlpha = std::atof(argv[2]);

    // First guess
    VectorXd initGuess = StdVecToEig(integrationParam.ODE.refParam);

    // Initialize MCMC
    MCMC algorithm(param, &priorObj, &likObj, observations,
                   nMCMC, RAM, initGuess, desiredAlpha);
    std::vector<VectorXd> result = algorithm.compute(&forwardGenerator);

    // Output
    std::ofstream output(std::string(DATA_PATH) + argv[3] + ".txt",
                         std::ofstream::out | std::ofstream::trunc);

    for (auto it : result) {
        output << it.transpose() << std::endl;
    }

    output.close();

    return 0;
}