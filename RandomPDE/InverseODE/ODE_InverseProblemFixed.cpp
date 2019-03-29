#include <MCMC.hpp>
#include <fstream>
#include <GetPot>
#include "ReadObservations.hpp"

int main(int argc, char* argv[])
{
    GetPot parser(argc, argv);
    Proposals proposal;

    double noise = 0.1,
            h = 0.1,
            p = 1.5,
            proposalStdDev = 0.1;
    int nObs = 10, nMCMC = 10000, nMC = 10;
    std::string outputFileName, obsFileName;
    bool prob = false, prob2 = false;

    if (parser.search("-noise"))
        noise = parser.next(noise);
    if (parser.search("-h"))
        h = parser.next(h);
    if (parser.search("-nObs"))
        nObs = parser.next(nObs);
    if (parser.search("-nMCMC"))
        nMCMC = parser.next(nMCMC);
    if (parser.search("-propVar"))
        proposalStdDev = parser.next(proposalStdDev);
    if (parser.search("-outputFile"))
        outputFileName = parser.next(" ");
    if (parser.search("-prob"))
        prob = true;
    if (parser.search("-obsFile"))
        obsFileName = parser.next(" ");
    if (parser.search("-p"))
        p = parser.next(p);
    if (parser.search("-prob2"))
        prob2 = true;
    if (parser.search("-nMC"))
        nMC = parser.next(nMC);

    // Initialize random seeds
    std::default_random_engine generatorOne{(unsigned int) time(NULL)};
    std::default_random_engine generatorTwo{(unsigned int) time(NULL) + 10};
    std::normal_distribution<double> noiseDist(0.0, noise);

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
    double T = tObs.back();


    // Method for the inverse solver
    Butcher tableau(EULERFORWARD, EXPLICIT, 0);

    // Initialize the MCMC machinery
    Posterior* posterior;
    proposal = Proposals(proposalStdDev);

    if (!prob && !prob2) {
        posterior = new RKPosteriorIC(h, observations, tObs, noise, ODE, tableau);
        MCMC mcmc(ODE.initialCond, &proposal, posterior,
                  static_cast<unsigned long>(nMCMC));
        std::vector<VectorXd> samples;

        VectorXd initGuess = VectorXd::Zero(ODE.size);
        samples = mcmc.compute(&generatorOne, &generatorTwo);
        std::ofstream output(DATA_PATH + outputFileName + ".txt", std::ofstream::out | std::ofstream::trunc);
        for(auto it : samples) {
            output << it.transpose() << std::endl;
        }
        output.close();
    } else if (prob && !prob2) {
        for (int i = 0; i < 100; i++) {
            posterior = new RKProbPosteriorFixed(h, T, ODE.initialCond, observations, tObs, noise, ODE, tableau, p, i);

            MCMC mcmc(ODE.refParam, &proposal, posterior, static_cast<unsigned long>(nMCMC));
            std::vector<VectorXd> samples;

            samples = mcmc.compute(&generatorOne, &generatorTwo);
            std::ofstream output(DATA_PATH + outputFileName + std::to_string(i) + ".txt", std::ofstream::out | std::ofstream::trunc);
            for (auto it : samples) {
                output << it.transpose() << std::endl;
            }
            output.close();
        }
    } else {
        posterior = new RKProbPosteriorIC(h, observations, tObs, noise, ODE, tableau, nMC, p);

        MCMC mcmc(ODE.initialCond, &proposal, posterior, static_cast<unsigned long>(nMCMC));
        std::vector<VectorXd> samples;

        samples = mcmc.compute(&generatorOne, &generatorTwo);
        std::ofstream output(DATA_PATH + outputFileName + ".txt", std::ofstream::out | std::ofstream::trunc);
        for (auto it : samples) {
            output << it.transpose() << std::endl;
        }
        output.close();
    }


    // delete posterior;

    return 0;
}