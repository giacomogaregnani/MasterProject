#include <RandomTimeStep.hpp>
#include <MCMC.hpp>
#include <UtilitiesGG.hpp>
#include <GetPot>

// GLOBAL VARIABLES

// Observations
VectorXd observations;
std::vector<double> obsTimes;
std::vector<unsigned long> obsIndex = {};
double obsNoise;

// Prior
VectorXd priorMean;

// Integration parameters
odeDef ODE;
Butcher tableau(IMPEULER, IMPLICIT, 0);
unsigned long N;
double h;

// RK Solver
RungeKuttaRandomH solver;
RungeKutta detSolver;
std::vector<std::shared_ptr<RungeKuttaRandomH>> MCSolver;

// Random number generator
std::default_random_engine generator{(unsigned int) time(NULL)};

// MCMC Parameters
int nMC;
int nMCMC;
int nThreads;

// POSTERIOR DEFINITION FOR MCMC



// =============================
int main(int argc, char* argv[])
{
    // Parse input parameters
    GetPot parser(argc, argv);
    std::string obsFileName = "observations", outFileName = "samples";
    bool noisy = false, RAM = false;
    double propVar = 0.1;
    nMC = 10;
    nMCMC = 10000;
    h = 0.1;

    if (parser.search("-obsFile"))
        obsFileName = parser.next("observations");
    if (parser.search("-nMC"))
        nMC = parser.next(nMC);
    if (parser.search("-nMCMC"))
        nMCMC = parser.next(nMCMC);
    if (parser.search("-outFile"))
        outFileName = parser.next("samples");
    if (parser.search("-h"))
        h = parser.next(h);
    if (parser.search("-noisy"))
        noisy = true;
    if (parser.search("-RAM"))
        RAM = true;
    if (parser.search("-prop"))
        propVar = parser.next(0.1);

    MCSolver.resize(nMC);
    if (nMC > 25)
        nThreads = 25;
    else
        nThreads = nMC;

    // Load observations
    std::string obsFile = std::string(DATA_PATH) + "HeatEquation/" + obsFileName + ".txt";
    loadObservationsHeat(observations, obsFile, &ODE.size, &obsNoise);

    // ODE
    ODE.ode = HEAT;
    setProblem(&ODE, ODE.size);
    std::vector<double> theta = {0.1};

    std::cout << observations.transpose() << std::endl;
    std::cout << obsNoise << std::endl;
    std::cout << ODE.size << std::endl;

    return 0;

    /* // Final time
    double T = 1;
    N = static_cast<unsigned long>(T / h);

    // Prior mean
    priorMean = VectorXd::Zero(ODE.refParam.size());

    // MCMC parameters
    VectorXd eigRefParam = StdVecToEig(ODE.refParam);
    std::normal_distribution<double>::param_type proposalParam(0.0, propVar);
    std::default_random_engine generatorMCMC{(unsigned int) time(NULL) + 1};

    // Probabilistic mcmc
    MCMC mcmc(eigRefParam, proposalParam, &posterior,
              (unsigned long) nMCMC, RAM, 0.234);
    std::vector<VectorXd> probResult = mcmc.compute(&generatorMCMC, noisy);

    // Deterministic mcmc
    MCMC detMcmc(eigRefParam, proposalParam, &detPosterior,
                 (unsigned long) nMCMC, true, 0.234);
    std::vector<VectorXd> detResult = detMcmc.compute(&generatorMCMC, false);

    // Write results on file
    std::ofstream detOutput(std::string(DATA_PATH) + "MCMC/" + outFileName + "Det.txt",
                            std::ofstream::out | std::ofstream::trunc);
    std::ofstream probOutput(std::string(DATA_PATH) + "MCMC/" + outFileName + "Prob.txt",
                             std::ofstream::out | std::ofstream::trunc);

    for (auto it : probResult) {
        probOutput << it.transpose() << std::endl;
    }

    for (auto it : detResult) {
        detOutput << it.transpose() << std::endl;
    }

    detOutput.close();
    probOutput.close();

    return 0; */
}
