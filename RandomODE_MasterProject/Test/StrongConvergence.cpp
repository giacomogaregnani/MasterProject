#include <RandomTimeStep.hpp>
#include <iomanip>
#include <GetPot>

int main(int argc, char* argv[])
{
    // Parse input parameters with GetPot
    GetPot parser(argc, argv);

    double h = 0.01;
    std::string outputName = "outputStrong";
    double p = 1.5;
    double T = 20;
    int nTrajectories = 1;
    int nExp = 1;
    int detOrder = 1;

    if (parser.search("-h"))
        h = parser.next(h);
    if (parser.search("-nExp"))
        nExp = parser.next(nExp);
    if (parser.search("-output"))
        outputName = parser.next("output");
    if (parser.search("-p"))
        p = parser.next(p);
    if (parser.search("-T"))
        T = parser.next(T);
    if (parser.search("-nTrajectories"))
        nTrajectories = parser.next(nTrajectories);
    if (parser.search("-detOrder"))
        detOrder = parser.next(detOrder);

    // Order 1, 2, 3, 4
    methods RKMethod = EULERFORWARD;
    switch (detOrder) {
        case 1:
            RKMethod = EULERFORWARD;
            break;
        case 2:
            RKMethod = EXPTRAPEZ;
            break;
        case 3:
            RKMethod = KU3;
        case 4:
            RKMethod = RK4;
    }

    // Which ODE ? Which RK method?
    odeDef ODE;
    ODE.ode = FITZNAG;
    setProblem(&ODE);
    int N = static_cast<int> (std::round(T / h));
    Butcher tableau(RKMethod, EXPLICIT, 0);

    // Initialize a random generator (use clock to make seed "random")
    std::default_random_engine generator{(unsigned int) time(NULL)};

    // Output file
    std::string outputFilename = std::string(DATA_PATH) + outputName + ".txt";
    std::ofstream output;
    output.open(outputFilename, std::ofstream::out | std::ofstream::trunc);

    // Reference solution
    double hRef = 1e-5;
    int NRef = static_cast<int> (std::round(T / hRef));
    RungeKutta refSolver(ODE, ODE.refParam, tableau);
    VectorXd refSolution = ODE.initialCond;
    for (int i = 0; i < NRef; i++) {
        refSolution = refSolver.oneStep(refSolution, hRef);
    }

    // Probabilistic solution

    // Time step loop
    for (int j = 0; j < nExp; j++) {
        std::cout << "Computing for h = " << h << std::endl << " =============== " << std::endl;

        RungeKuttaRandomH probSolver(&generator, ODE, ODE.refParam, tableau, h, p);

        double err = 0;
        // Monte Carlo loop
        for (int k = 0; k < nTrajectories; k++) {
            VectorXd probSolution = ODE.initialCond;

            // Time integration loop
            for (int i = 0; i < N; i++)
                probSolution = probSolver.oneStep(probSolution);

            err += (refSolution - probSolution).norm();
        }
        err /= nTrajectories;

        output << h << "\t" << err << std::endl;
        h /= 2; N *= 2;
    }

    output.close();
    return 0;
}