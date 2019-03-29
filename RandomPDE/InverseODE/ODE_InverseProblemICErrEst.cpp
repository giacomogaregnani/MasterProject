#include <MCMC.hpp>
#include <fstream>
#include <GetPot>
#include <iomanip>
#include "ReadObservations.hpp"

int main(int argc, char* argv[]) {
    GetPot parser(argc, argv);
    Proposals proposal;

    double  noise = 0.1,
            h = 0.1,
            proposalStdDev = 0.1;
    int nObs = 1, nMCMC = 10000, nErrEst = 10, nLoop = 5;
    std::string outputFileName, obsFileName;
    bool isGauss = false;

    if (parser.search("-noise"))
        noise = parser.next(noise);
    if (parser.search("-h"))
        h = parser.next(h);
    if (parser.search("-nMCMC"))
        nMCMC = parser.next(nMCMC);
    if (parser.search("-propVar"))
        proposalStdDev = parser.next(proposalStdDev);
    if (parser.search("-outputFile"))
        outputFileName = parser.next(" ");
    if (parser.search("-obsFile"))
        obsFileName = parser.next(" ");
    if (parser.search("-nObs"))
        nObs = parser.next(nObs);
    if (parser.search("-Gauss"))
        isGauss = true;
    if (parser.search("-nErr"))
        nErrEst = parser.next(nErrEst);
    if (parser.search("-L"))
        nLoop = parser.next(nLoop);

    // Initialize random seeds
    std::default_random_engine generatorOne{(unsigned int) time(NULL)};
    std::default_random_engine generatorTwo{(unsigned int) time(NULL) + 10};

    // Model ODE
    odeDef ODE;
    ODE.ode = TEST1D;
    setProblem(&ODE);

    // Import observations
    obsFileName = DATA_PATH + obsFileName + ".txt";
    std::vector<VectorXd> observations = {};
    std::vector<double> tObs = {};
    ReadObservations(tObs, observations, (unsigned int) nObs,
                     (unsigned int) ODE.size, obsFileName);

    // Method for the inverse solver
    Butcher tableau;
    if (isGauss) {
        tableau = Butcher(SYMPEULER, SEPAR);
    } else {
        tableau = Butcher(EULERFORWARD, EXPLICIT);
    }

    double hRef = 1e-3;
    auto NRef = static_cast<int>(tObs.back() / hRef);
    auto N = static_cast<int>(tObs.back() / h);
    Butcher refTableau(RK4, EXPLICIT);
    RungeKutta refSolver(ODE, refTableau);
    RungeKutta detSolver(ODE, tableau);
    std::normal_distribution<double> prior(0.0, 1.0);
    std::vector<VectorXd> samples = {VectorXd::Zero(ODE.size)};

    for (int k = 0; k < nLoop; k++) {
        std::uniform_int_distribution<int> unifDist(0, static_cast<int>(samples.size())-1);
        std::vector<VectorXd> errors(nErrEst, VectorXd::Zero(ODE.size));

        for (int i = 0; i < nErrEst; i++) {
            VectorXd newIC(ODE.size);

            if (k == 0) {
                for (int j = 0; j < ODE.size; j++) {
                    newIC(j) = prior(generatorOne);
                }
            } else {
                int idx = unifDist(generatorOne);
                newIC = samples[idx];
            }

            // Solve  the equation for this new initial condition
            VectorXd fineSolution = newIC;
            for (int j = 0; j < NRef; j++) {
                fineSolution = refSolver.oneStep(hRef, fineSolution, ODE.refParam);
            }
            VectorXd coarseSolution = newIC;
            for (int j = 0; j < N; j++) {
                coarseSolution = detSolver.oneStep(h, coarseSolution, ODE.refParam);
            }

            errors[i] = fineSolution - coarseSolution;
        }



        // Initialize the MCMC machinery
        std::normal_distribution<double> noiseDist(0.0, noise);
        proposal = Proposals(proposalStdDev);
        RKEEPosteriorIC posterior(h, observations, errors, tObs, noise, ODE, tableau);

        // Compute the posterior with MCMC
        VectorXd initGuess = samples.back();

        int effNMCMC = nMCMC / nLoop * (k < nLoop - 1) + nMCMC * (k == nLoop - 1);
        MCMC mcmc(initGuess, &proposal, &posterior, effNMCMC);
        samples = mcmc.compute(&generatorOne, &generatorTwo);
    }

    std::ofstream output(DATA_PATH + outputFileName + ".txt", std::ofstream::out | std::ofstream::trunc);
    for (auto it : samples) {
        output << std::fixed << std::setprecision(20) << it.transpose() << std::endl;
    }
    output.close();

    return 0;
}