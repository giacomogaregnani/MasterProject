#include <RungeKuttaSolver.hpp>
#include <GetPot>
#include <random>

int main(int argc, char* argv[])
{
    GetPot parser(argc, argv);

    double noise = 0.1,
           T = 1,
           h = 1e-3;
    int nObs = 10;
    std::string outputFileName;

    if (parser.search("-h"))
        h = parser.next(h);
    if (parser.search("-noise"))
        noise = parser.next(noise);
    if (parser.search("-T"))
        T = parser.next(T);
    if (parser.search("-nObs"))
        nObs = parser.next(nObs);
    if (parser.search("-outputFile"))
        outputFileName = parser.next(" ");

    // Initialize random seeds
    std::default_random_engine generatorOne{2018};
    std::normal_distribution<double> noiseDist(0.0, noise);

    // Model ODE
    odeDef ODE;
    ODE.ode = HENHEIL;
    setProblem(&ODE);

    // Generate observations
    std::vector<VectorXd> observations = {};
    std::vector<double> tObs = {};
    VectorXd refTheta = ODE.refParam;

    double tIntObs = T / nObs;
    for (int i = 1; i < nObs + 1; i++)
        tObs.push_back(i * tIntObs);
    Butcher refTableau(EULERFORWARD, EXPLICIT);
    RungeKutta refSolver(ODE, refTableau);

    VectorXd refSolution = ODE.initialCond;
    VectorXd vecNoise(ODE.size);
    tObs.insert(tObs.begin(), 0.0);

    for (unsigned int i = 0; i < tObs.size() - 1; i++) {
        double dT = tObs[i+1] - tObs[i];
        auto nSteps = static_cast<unsigned long>(std::round(dT / h));

        for (unsigned long j = 0; j < nSteps; j++) {
            refSolution = refSolver.oneStep(h, refSolution, refTheta);
        }

        for (int j = 0; j < ODE.size; j++) {
            vecNoise(j) = noiseDist(generatorOne);
        }

        observations.push_back(refSolution + vecNoise);
        vecNoise = VectorXd::Zero(ODE.size);
    }
    tObs.erase(tObs.begin());

    std::ofstream output(DATA_PATH + outputFileName + ".txt", std::ofstream::out | std::ofstream::trunc);
    for (unsigned int i = 0; i < tObs.size(); i++) {
        output << tObs[i] << "\t" << observations[i].transpose() << std::endl;
    }
    output.close();


    return 0;
}