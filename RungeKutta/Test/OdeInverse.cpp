#include <RandomTimeStep.hpp>
#include <MCMC.hpp>
#include <UtilitiesGG.hpp>
#include <GetPot>

// GLOBAL VARIABLES

// Observations
std::vector<VectorXd> observations;
std::vector<double> obsTimes;
std::vector<unsigned long> obsIndex = {};
double obsNoise;

// Prior
VectorXd priorMean;

// Integration parameters
odeDef ODE;
Butcher tableau(EULERFORWARD, EXPLICIT, 0);
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

double MClikelihood(VectorXd& theta)
{
    /*VectorXd expTheta(theta.size());
    for (int i = 0; i < theta.size(); i++) {
        expTheta(i) = exp(theta(i));
    }*/

    std::vector<double> stdTheta = EigVecToStdVec(theta);
    std::vector<double> likMC(nMC, 0.0);
    size_t k;
    int j;

    for (j = 0; j < nMC; j++) {
        MCSolver[j] = std::make_shared<RungeKuttaRandomH>(&generator, ODE, stdTheta, tableau, h, 1.5);
    }

    #pragma omp parallel for num_threads(nThreads) private(j, k)
    for (j = 0; j < nMC; j++) {

        VectorXd solution = ODE.initialCond;
        k = 0;

        for (unsigned long i = 0; i < N; i++) {
            solution = MCSolver[j]->oneStep(solution);
            if (i + 1 == obsIndex[k]) {
                likMC[j] += -0.5 / (obsNoise * obsNoise) * (solution - observations[k]).dot(solution - observations[k]);
                ++k;
            }
        }
    }

    std::vector<double>::iterator maxLikIt = std::max_element(likMC.begin(), likMC.end());
    double maxLik = *maxLikIt;
    likMC.erase(maxLikIt, maxLikIt);

    double sum = 0;
    for (auto it : likMC) {
        sum += exp(it - maxLik);
    }

    return maxLik + std::log(1 + sum) - std::log(static_cast<double>(nMC));
}

double detLikelihood(VectorXd& theta)
{
    /* VectorXd expTheta(theta.size());
    for (int i = 0; i < theta.size(); i++) {
        expTheta(i) = exp(theta(i));
    } */

    std::vector<double> stdTheta = EigVecToStdVec(theta);
    detSolver = RungeKutta(ODE, stdTheta, tableau);

    VectorXd solution = ODE.initialCond;

    double lik = 0;

    unsigned int k = 0;
    for (unsigned long i = 0; i < N; i++) {
        solution = detSolver.oneStep(solution, h);
        if (i + 1 == obsIndex[k]) {
            lik += -0.5 / (obsNoise * obsNoise) * (solution - observations[k]).dot(solution - observations[k]);
            ++k;
        }
    }

    return lik;
}

double prior(VectorXd& theta)
{
    return -0.5 * (theta - priorMean).dot(theta - priorMean);
}

double posterior(VectorXd& theta)
{
    return MClikelihood(theta) + prior(theta);
}

double detPosterior(VectorXd& theta)
{
    return detLikelihood(theta) + prior(theta);
}

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

    // ODE
    ODE.ode = LORENZ;
    setProblem(&ODE);

    // Load observations
    std::string obsFile = std::string(DATA_PATH) + "MCMC/" + obsFileName + ".txt";
    loadObservations(obsTimes, observations, obsFile, ODE.size, &obsNoise);

    // Final time
    double T = obsTimes.back();
    N = static_cast<unsigned long>(T / h);

    // Compute observation steps
    for (auto it : obsTimes) {
        obsIndex.push_back(static_cast<unsigned long> (it / h));
    }

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

    return 0;
}
