#include <RandomTimeStep.hpp>
#include <iomanip>
#include <GetPot>

int main(int argc, char* argv[])
{
    // Parse input parameters with GetPot
    GetPot parser(argc, argv);

    double h = 0.01;
    std::string outputName = "output";
    double p = 1.5;
    double T = 20;
    int nTrajectories = 1;

    if (parser.search("-h"))
        h = parser.next(h);
    if (parser.search("-output"))
        outputName = parser.next("output");
    if (parser.search("-p"))
        p = parser.next(p);
    if (parser.search("-T"))
        T = parser.next(T);
    if (parser.search("-nTrajectories"))
        nTrajectories = parser.next(nTrajectories);

    // Which ODE ? Which RK method?
    odeDef ODE;
    ODE.ode = FITZNAG;
    setProblem(&ODE);
    int N = static_cast<int> (std::round(T / h));
    Butcher tableau(EULERFORWARD, EXPLICIT, 0);

    // Initialize a random generator (use clock to make seed "random")
    std::default_random_engine generator{(unsigned int) time(NULL)};

    // Output file
    std::string outputFilename = std::string(DATA_PATH) + outputName + ".txt";
    std::ofstream output;
    output.open(outputFilename, std::ofstream::out | std::ofstream::trunc);

    // Deterministic solution
    RungeKutta detSolver(ODE, ODE.refParam, tableau);
    VectorXd detSolution = ODE.initialCond;
    output << ODE.initialCond.transpose() << std::endl;
    for (int i = 0; i < N; i++) {
        detSolution = detSolver.oneStep(detSolution, h);
        output << detSolution.transpose() << std::endl;
    }

    // Probabilistic solution
    RungeKuttaRandomH probSolver(&generator, ODE, ODE.refParam, tableau, h, p);
    for (int k = 0; k < nTrajectories; k++) {
        VectorXd probSolution = ODE.initialCond;

        // Time integration loop
        output << ODE.initialCond.transpose() << std::endl;
        for (int i = 0; i < N; i++) {
            probSolution = probSolver.oneStep(probSolution);
            output << probSolution.transpose() << std::endl;
        }
    }

    output.close();
    return 0;
}