#include <RandomTimeStep.hpp>
#include <fstream>
#include <iomanip>
#include <GetPot>

double phi(VectorXd& v)
{
    return v.dot(v);
}

int main(int argc, char* argv[])
{
    GetPot parser(argc, argv);

    double h = 0.1,
            T = 50;
    int nMC = 10, nExp = 1;

    std::string outputFile;

    if (parser.search("-h"))
        h = parser.next(h);
    if (parser.search("-T"))
        T = parser.next(T);
    if (parser.search("-output"))
        outputFile = parser.next(" ");
    if (parser.search("-nMC"))
        nMC = parser.next(nMC);
    if (parser.search("-nExp"))
        nExp = parser.next(nExp);

    // ODE
    odeDef ODE;
    ODE.ode = LORENZ;
    setProblem(&ODE);

    // Numerical method
    Butcher tableau(EULERFORWARD, EXPLICIT, 0);

    // Integration parameters
    auto N = static_cast<unsigned int>(T / h);
    double hDet = h / 100;
    double TDet = T * 10;
    auto NDet = static_cast<unsigned int>(TDet / hDet);

    // Random seed
    std::default_random_engine generator{(unsigned int) time(NULL)};

    std::string outputFileName = std::string(DATA_PATH) + outputFile + ".txt";
    std::ofstream output(outputFileName, std::ofstream::out | std::ofstream::trunc);

    // Deterministic (Longer Time)
    RungeKutta detMethod(ODE, tableau);
    VectorXd solution = ODE.initialCond;
    double phiDet = 0;
    for (unsigned int i = 0; i < NDet; i++) {
        solution = detMethod.oneStep(hDet, solution, ODE.refParam);
        phiDet += phi(solution);
    }
    phiDet /= NDet;
    output << phiDet << std::endl;

    // Initialization
    RungeKuttaRandomH Method(&generator, ODE, tableau, h, 1.5);
    double phiProb;
    double phiProbMean;

    for (unsigned int k = 0; k < nExp; k++) {
        phiProbMean = 0.0;
        for (unsigned int l = 0; l < 50; l++) {
            phiProb = 0.0;
            for (unsigned int j = 0; j < nMC; j++) {
                solution = ODE.initialCond;
                for (unsigned int i = 0; i < N; i++) {
                    solution = Method.oneStep(solution, ODE.refParam);
                }
                phiProb += phi(solution);
            }
            phiProbMean += phiProb / nMC;
        }
        phiProbMean /= 50;
        output << nMC << "\t" << phiProbMean << std::endl;

        nMC *= 2;
    }

    output.close();

    return 0;
}
