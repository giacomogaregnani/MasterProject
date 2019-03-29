#include <MCMC.hpp>
#include <fstream>
#include <GetPot>

double f(double x)
{
    // return exp(-100*std::pow(x - 1/2, 2)) * sin(12*PI*x)* std::pow(200*x - 100, 2) - 200*exp(-100*std::pow(x - 1/2, 2))*sin(12*PI*x) -
    //       144*PI*PI*exp(-100*std::pow(x - 1/2, 2))*sin(12*PI*x) - 24*PI*exp(-100*std::pow(x - 1/2, 2))*cos(12*PI*x)*(200*x - 100);
    return sin(2.0 * PI * x);
}

VectorXd buildField(VectorXd& param, double h)
{
    auto N = static_cast<int>(std::round(1.0 / h));
    VectorXd field = VectorXd::Zero(N+1);

    for (int i = 0; i < N+1; i++) {
        if (h*i <= 0.4 && h*i >= 0.2) {
            field(i) = std::log(1.0 + param(0));
        } else if (h*i <= 0.8 && h*i >= 0.6) {
            field(i) = std::log(1.0 + param(1));
        } else {
            field(i) = 0.0;
        }
    }

    for (int i = 0; i < N+1; i++)
        field(i) = std::exp(field(i));

    return field;
}

int main(int argc, char* argv[])
{
    GetPot parser(argc, argv);
    KarhunenLoeve KL;
    Proposals proposal;

    double  noise = 0.1,
            h = 0.1,
            proposalStdDev = 0.1;
    int     nObs = 9,
            nMCMC = 10000,
            nChains = 1;
    bool    saveSample = false,
            obsGiven = false;
    std::string outputFileName,
            obsFileName;

    if (parser.search("-noise"))
        noise = parser.next(noise);
    if (parser.search("-h"))
        h = parser.next(h);
    if (parser.search("-nMCMC"))
        nMCMC = parser.next(nMCMC);
    if (parser.search("-propVar"))
        proposalStdDev = parser.next(proposalStdDev);
    if (parser.search("-nObs"))
        nObs = parser.next(nObs);
    if (parser.search("-outputFile"))
        outputFileName = parser.next(" ");
    if (parser.search("-nChains"))
        nChains = parser.next(nChains);
    if (parser.search("-saveSample"))
        saveSample = true;
    if (parser.search("-obsGiven"))
        obsGiven = true;

    auto N = static_cast<int>(1.0 / h);
    VectorXd meanField = VectorXd::Zero(N+1);

    // Initialize the proposal structure and the prior parameters
    proposal = Proposals(proposalStdDev);

    // Initialize random seeds
    std::default_random_engine generatorOne{(unsigned int) time(nullptr)};
    std::default_random_engine generatorTwo{(unsigned int) time(nullptr) + 10};

    // Generate observations from true value and a reference solution
    std::fstream obsFile;
    VectorXd uRef;
    VectorXd xObs;
    xObs.setLinSpaced(nObs+2, 0.0, 1.0);
    xObs = xObs.segment(1, nObs);
    VectorXd observations(nObs);
    VectorXd trueParam(2);
    trueParam << 0.3, -0.3;

    double hRef = 1e-3;
    Solver refSolver(f, 0, 1, hRef, 0, 0);
    auto NRef = static_cast<int>(std::round(1 / hRef));

    // Generate the true value of the conductivity
    if (obsGiven) {
        std::ifstream obsInput(DATA_PATH + std::string("obsfin.txt"), std::ifstream::in);
        for (int i = 0; i < nObs; i++)
            obsInput >> observations(i);
        obsInput.close();
    }
    else {
        std::ofstream obsOutput(DATA_PATH + std::string("obsfin.txt"), std::ofstream::out | std::ofstream::trunc);
        std::cout << "Generating true conductivity ..." << std::endl;
        VectorXd kappaRef(NRef + 1);
        kappaRef = buildField(trueParam, hRef);

        // Compute the observations
        std::cout << "Computing observations ..." << std::endl;
        std::normal_distribution<double> noiseDist(0.0, noise);
        uRef = refSolver.solve(kappaRef);

        observations = refSolver.evaluate(xObs);
        for (int i = 0; i < nObs; i++) {
            observations(i) += noiseDist(generatorOne);
            obsOutput << observations(i) << std::endl;
        }

    }

    // DETERMINISTIC COMPUTATION
    // Initialize the posterior machinery
    std::vector<double> BC = {0.0, 0.0};
    FEMPosteriorFin detPosterior(&f, h, BC, &buildField, observations, xObs, noise);

    // Compute the posterior with MCMC
    VectorXd initGuess = trueParam; // VectorXd::Zero(trueParam.size());
    MCMC mcmc(initGuess, &proposal, &detPosterior, static_cast<unsigned long>(nMCMC));
    std::vector<VectorXd> samples;
    std::cout << "Computing sample ..." << std::endl;
    samples = mcmc.compute(&generatorOne, &generatorTwo, false, false);

    // Write on file the raw samples
    if (saveSample) {
        std::ofstream rawOutput(DATA_PATH + outputFileName + "coeffs.txt", std::ofstream::out | std::ofstream::trunc);
        for (auto it : samples)
            rawOutput << it.transpose() << std::endl;
        rawOutput.close();
    }

    // PROBABILISTIC COMPUTATION
    // Initialize the posterior machinery
    FEMProbPosteriorFin probPosterior(&f, h, BC, 1.0, &buildField, observations, xObs, noise);
    MCMC mcmcProb(initGuess, &proposal, &probPosterior, static_cast<unsigned long>(nMCMC));

    std::vector<VectorXd> samplesProb, samplesTemp;
    for (int i = 0; i < nChains; i++) {
        mcmcProb.eraseSample();
        samplesTemp = mcmcProb.compute(&generatorOne, &generatorTwo, false, false);
        samplesProb.insert(samplesProb.end(), samplesTemp.begin(), samplesTemp.end());
        probPosterior.resetMesh();
    }

    // Write on file the raw samples
    if (saveSample) {
        std::ofstream rawProbOutput(DATA_PATH + outputFileName + "coeffsProb.txt",
                                    std::ofstream::out | std::ofstream::trunc);
        for (auto it : samplesProb)
            rawProbOutput << it.transpose() << std::endl;
        rawProbOutput.close();
    }

    return 0;
}