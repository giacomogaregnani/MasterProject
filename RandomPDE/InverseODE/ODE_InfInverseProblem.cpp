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
            T = 1,
            p = 1.5,
            proposalStdDev = 0.1;
    int nObs = 10,
            nMC = 10,
            nMCMC = 10000,
            nRKC = 10,
            nKL = 10;
    std::string outputFileName, obsFileName;
    bool prob = false, noisy = false;

    if (parser.search("-noise"))
        noise = parser.next(noise);
    if (parser.search("-h"))
        h = parser.next(h);
    if (parser.search("-T"))
        T = parser.next(T);
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
    if (parser.search("-nMC"))
        nMC = parser.next(nMC);
    if (parser.search("-obsFile"))
        obsFileName = parser.next(" ");
    if (parser.search("-noisy"))
        noisy = true;
    if (parser.search("-p"))
        p = parser.next(p);
    if (parser.search("-nKL"))
        nKL = parser.next(nKL);
    if (parser.search("-nRKC"))
        nRKC = parser.next(nRKC);

    // Initialize random seeds
    std::default_random_engine generatorOne{(unsigned int) time(NULL)};
    std::default_random_engine generatorTwo{(unsigned int) time(NULL) + 10};
    std::normal_distribution<double> noiseDist(0.0, noise);

    // Model ODE
    odeDef ODE;
    ODE.ode = BRUSS;
    setProblem(&ODE);

    // Import observations
    obsFileName = DATA_PATH + obsFileName + ".txt";
    std::vector<VectorXd> observations = {};
    std::vector<double> tObs = {};
    ReadObservations(tObs, observations, (unsigned int) nObs,
                     (unsigned int) ODE.size, obsFileName);


    // Method for the inverse solver
    Butcher tableau(RKC, STABEXP, nRKC);

    // Initialize the MCMC machinery
    std::vector<double> proposalParam(1, proposalStdDev);
    proposal = Proposals(proposalStdDev, proposalParam);

    Posterior* posterior;
    VectorXd KLMean = VectorXd::Ones(ODE.size / 2 + 2);

    if (!prob) {
        posterior = new RKPosteriorIC(h, T, observations, tObs, noise, ODE, tableau, KLMean);
    } else {
        posterior = new RKProbPosteriorIC(h, T, observations, tObs, noise, ODE, tableau, KLMean, nMC, p);
    }

    // Compute the posterior with MCMC
    VectorXd initGuess = VectorXd::Zero(nKL);
    MCMC mcmc(initGuess, &proposal, posterior, static_cast<unsigned long>(nMCMC));
    std::vector<VectorXd> samples;
    samples = mcmc.compute(&generatorOne, &generatorTwo, noisy);

    std::ofstream outputTheta(DATA_PATH + outputFileName + "theta.txt", std::ofstream::out | std::ofstream::trunc);
    for(auto it : samples) {
        outputTheta << it.transpose() << std::endl;
    }
    outputTheta.close();

    std::ofstream output(DATA_PATH + outputFileName + ".txt", std::ofstream::out | std::ofstream::trunc);
    for(auto it : samples) {
        output << (*posterior).coeffToField(it).transpose() << std::endl;
    }
    output.close();

    return 0;
}
