#include <MCMC.hpp>
#include <fstream>
#include <GetPot>
#include <iomanip>
#include "ReadObservations.hpp"

int main(int argc, char* argv[])
{
    GetPot parser(argc, argv);
    Proposals proposal;

    double noise = 0.1,
           h = 0.1,
           p = 1.5,
           proposalStdDev = 0.1;
    int    nObs = 1,
           nMC = 10,
           nMCMC = 10000;
    std::string outputFileName, obsFileName;
    bool prob = false, noisy = false, isGauss = false;

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
    if (parser.search("-prob"))
        prob = true;
    if (parser.search("-nMC"))
        nMC = parser.next(nMC);
    if (parser.search("-obsFile"))
        obsFileName = parser.next(" ");
    if (parser.search("-nObs"))
        nObs = parser.next(nObs);
    if (parser.search("-noisy"))
        noisy = true;
    if (parser.search("-p"))
        p = parser.next(p);
    if (parser.search("-Gauss"))
        isGauss = true;

    // Initialize random seeds
    std::default_random_engine generatorOne{(unsigned int) time(NULL)};
    std::default_random_engine generatorTwo{(unsigned int) time(NULL) + 10};

    // Model ODE
    odeDef ODE;
    ODE.ode = HENHEIL;
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
        tableau = Butcher(STORMVER, SEPAR);
    } else {
        tableau = Butcher(EULERFORWARD, EXPLICIT);
    }

    std::normal_distribution<double> noiseDist(0.0, noise);

    // Initialize the MCMC machinery
    Posterior* posterior;
    proposal = Proposals(proposalStdDev);

    if (!prob) {
        posterior = new RKPosteriorIC(h, observations, tObs, noise, ODE, tableau);
    } else {
        posterior = new RKProbPosteriorIC(h, observations, tObs, noise, ODE, tableau, nMC, p);
    }

    // Compute the posterior with MCMC
    // VectorXd initGuess = VectorXd::Zero(ODE.size);
    VectorXd initGuess = ODE.initialCond; //+ 0.1 * VectorXd::Random(ODE.size);
    MCMC mcmc(initGuess, &proposal, posterior,
              static_cast<unsigned long>(nMCMC));
    std::vector<VectorXd> samples;
    samples = mcmc.compute(&generatorOne, &generatorTwo, noisy);

    std::ofstream output(DATA_PATH + outputFileName + ".txt", std::ofstream::out | std::ofstream::trunc);
    for(auto it : samples) {
        output << std::fixed << std::setprecision(20) << it.transpose() << std::endl;
    }
    output.close();

    return 0;
}