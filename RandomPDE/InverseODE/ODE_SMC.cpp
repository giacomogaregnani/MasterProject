#include <SMC.hpp>
#include <SMCDet.hpp>
#include <fstream>
#include <GetPot>
#include "ReadObservations.hpp"

int main(int argc, char* argv[])
{
    GetPot parser(argc, argv);

    double noise = 0.1,
           h = 0.1,
           T = 1;
    int nObs = 10,
        nParticles = 10000;
    bool prob = false;
    std::string outputFileName, obsFileName;

    if (parser.search("-noise"))
        noise = parser.next(noise);
    if (parser.search("-h"))
        h = parser.next(h);
    if (parser.search("-T"))
        T = parser.next(T);
    if (parser.search("-nObs"))
        nObs = parser.next(nObs);
    if (parser.search("-nPart"))
        nParticles = parser.next(nParticles);
    if (parser.search("-outputFile"))
        outputFileName = parser.next(" ");
    if (parser.search("-obsFile"))
        obsFileName = parser.next(" ");
    if (parser.search("-prob"))
        prob = true;

    // Initialize random seeds
    std::default_random_engine generator{(unsigned int) time(NULL)};
    std::normal_distribution<double> noiseDist(0.0, noise);

    // Model ODE
    odeDef ODE;
    ODE.ode = FITZNAG;
    setProblem(&ODE);

    // Import observations
    obsFileName = DATA_PATH + obsFileName + ".txt";
    std::vector<VectorXd> observations = {};
    std::vector<double> tObs = {};
    ReadObservations(tObs, observations, (unsigned int) nObs,
                     (unsigned int) ODE.size, obsFileName);

    // Method for the inverse solver
    Butcher tableau(EULERFORWARD, EXPLICIT, 0);
    Butcher tableau2(EXPTRAPEZ, EXPLICIT, 0);
    double p = 1.5;

    // Compute the posterior with SMC
    VectorXd priorMean = VectorXd::Zero(ODE.refParam.size());
    SMC smc((unsigned long) nParticles, 0.97, ODE.initialCond, ODE, tableau, h, p, &generator, priorMean, 1.0);
    SMCDet smc2((unsigned long) nParticles, 0.97, ODE.initialCond, ODE, tableau, tableau2, h, &generator, priorMean, 1.0);

    std::vector<VectorXd> samples = {}, thetas = {};
    std::vector<double> weights = {};
    if (prob) {
        smc.compute(samples, thetas, weights, tObs, observations, noise);
    } else {
        smc2.compute(samples, thetas, weights, tObs, observations, noise);
    }

    std::ofstream output(DATA_PATH + outputFileName + ".txt", std::ofstream::out | std::ofstream::trunc);
    for (int i = 0; i < nParticles; i++) {
        output << weights[i] << " " << thetas[i].transpose() << std::endl;
    }
    output.close();

    return 0;
}