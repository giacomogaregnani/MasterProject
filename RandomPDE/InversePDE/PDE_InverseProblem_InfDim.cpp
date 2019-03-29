#include <MCMC.hpp>
#include <fstream>
#include <GetPot>

double f(double x)
{
    // return exp(-100*std::pow(x - 1/2, 2)) * sin(12*PI*x)* std::pow(200*x - 100, 2) - 200*exp(-100*std::pow(x - 1/2, 2))*sin(12*PI*x) -
    //       144*PI*PI*exp(-100*std::pow(x - 1/2, 2))*sin(12*PI*x) - 24*PI*exp(-100*std::pow(x - 1/2, 2))*cos(12*PI*x)*(200*x - 100);
    return sin(2.0 * PI * x);
}

int main(int argc, char* argv[])
{
    GetPot parser(argc, argv);
    KarhunenLoeve KL;
    Proposals proposal;

    double  noise = 0.1,
            h = 0.1,
            proposalStdDev = 0.1,
            gammaExpPrior = 0.1,
            deltaExpPrior = 0.1;
    int     nObs = 9,
            nMCMC = 10000,
            nEig = 2,
            nChains = 1;
    bool    disc = false,
            saveSample = false;
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
    if (parser.search("-gammaExp"))
        gammaExpPrior = parser.next(gammaExpPrior);
    if (parser.search("-deltaExp"))
        deltaExpPrior = parser.next(deltaExpPrior);
    if (parser.search("-nObs"))
        nObs = parser.next(nObs);
    if (parser.search("-outputFile"))
        outputFileName = parser.next(" ");
    if (parser.search("-nEig"))
        nEig = parser.next(nEig);
    if (parser.search("-nChains"))
        nChains = parser.next(nChains);
    if (parser.search("-discontinuous"))
        disc = true;
    if (parser.search("-saveSample"))
        saveSample = true;

    auto N = static_cast<int>(1.0 / h);
    VectorXd meanField = VectorXd::Zero(N+1);

    // Initialize the proposal structure and the prior parameters
    proposal = Proposals(proposalStdDev);
    std::vector<double> priorParam = {gammaExpPrior, deltaExpPrior};

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
    VectorXd trueParam(4);
    trueParam << 1.0, 1.0, 0.0, 1.0;

    double hRef = 1e-3;
    Solver refSolver(f, 0, 1, hRef, 0, 0);
    auto NRef = static_cast<int>(std::round(1 / hRef));

    // Generate the true value of the conductivity
    std::cout << "Generating true conductivity ..." << std::endl;
    VectorXd kappaRef(NRef+1);
    covKernelFcts covPrior = SQDEXP;
    if (disc) {
        for (int i = 0; i < NRef+1; i++) {
            if (hRef*i <= 0.4 && hRef*i >= 0.2) {
                kappaRef(i) = std::log(1.5);
            }
            if (hRef*i <= 0.8 && hRef*i >= 0.6) {
                kappaRef(i) = std::log(0.5);
            }
        }
    } else {
        VectorXd fieldMean = VectorXd::Zero(NRef + 1);
        KarhunenLoeve KarLoe(fieldMean, covPrior, NRef + 1, priorParam);
        kappaRef = KarLoe.KL(trueParam);
    }
    for (int i = 0; i < kappaRef.size(); i++)
        kappaRef(i) = std::exp(kappaRef(i));

    // Compute the observations
    std::cout << "Computing observations ..." << std::endl;
    std::normal_distribution<double> noiseDist(0.0, noise);
    uRef = refSolver.solve(kappaRef);

    observations = refSolver.evaluate(xObs);
    for (int i = 0; i < nObs; i++) {
        observations(i) += noiseDist(generatorOne);
    }

    // Print on a file the true conductivity
    std::ofstream outputTruth(DATA_PATH + outputFileName + "truth.txt", std::ofstream::out | std::ofstream::trunc);
    outputTruth << kappaRef.transpose() << std::endl;
    outputTruth.close();

    // DETERMINISTIC COMPUTATION
    // Initialize the posterior machinery
    std::vector<double> BC = {0.0, 0.0};
    FEMPosterior detPosterior(&f, h, BC, meanField, covPrior, priorParam, observations, xObs, noise);

    // Compute the posterior with MCMC
    if (nEig > N - 1)
        nEig = N - 1;
    VectorXd initGuess = VectorXd::Zero(nEig);
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

    // Compute pushed sample mean and output it on file
    std::ofstream output(DATA_PATH + outputFileName + ".txt", std::ofstream::out | std::ofstream::trunc);
    VectorXd meanResult = VectorXd::Zero(N+1);
    for (auto it : samples)
        meanResult += detPosterior.coeffToField(it);
    meanResult /= samples.size();
    output << meanResult.transpose() << std::endl;

    // Compute pushed sample variance
    MatrixXd MCMCVar = MatrixXd::Zero(N+1, N+1);
    VectorXd diff;
    VectorXd tmp;
    for (auto it : samples) {
        tmp = detPosterior.coeffToField(it);
        diff = meanResult - tmp;
        MCMCVar += diff * diff.transpose();
    }
    MCMCVar /= (samples.size() - 1);
    output << MCMCVar << std::endl;
    output.close();

    // PROBABILISTIC COMPUTATION
    // Initialize the posterior machinery
    FEMProbPosteriorMulti probPosterior(&f, h, BC, 1.0, meanField, covPrior, priorParam, observations, xObs, noise);
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

    // Compute sample mean and outputProb it on file
    std::ofstream outputProb(DATA_PATH + outputFileName + "prob.txt", std::ofstream::out | std::ofstream::trunc);
    VectorXd meanResultProb = VectorXd::Zero(N+1);
    for (auto it : samplesProb)
        meanResultProb += probPosterior.coeffToField(it);
    meanResultProb /= samplesProb.size();
    outputProb << meanResultProb.transpose() << std::endl;

    // Compute sample variance
    MatrixXd MCMCVarProb = MatrixXd::Zero(N+1, N+1);
    MatrixXd MCMCVarSampleProb = VectorXd::Zero(N+1);
    for (auto it : samplesProb) {
        tmp = probPosterior.coeffToField(it);
        diff = meanResultProb - tmp;
        MCMCVarProb += diff * diff.transpose();
    }
    MCMCVarProb /= (samplesProb.size() - 1);
    outputProb << MCMCVarProb << std::endl;
    outputProb.close();

    return 0;
}