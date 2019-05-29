#include <RandomTimeStep.hpp>
#include <fstream>
#include <iomanip>
#include <GetPot>
#include <iostream>

double phi(VectorXd& x)
{
    return x.dot(x);
}

int main(int argc, char* argv[])
{
    GetPot parser(argc, argv);

    double h = 0.1,
            T = 1,
            p = 1.0;
    int nExp = 1, nMC = 10;

    std::string outputFileName;

    if (parser.search("-h"))
        h = parser.next(h);
    if (parser.search("-T"))
        T = parser.next(T);
    if (parser.search("-output"))
        outputFileName = parser.next(" ");
    if (parser.search("-nExp"))
        nExp = parser.next(nExp);
    if (parser.search("-p"))
        p = parser.next(p);
    if (parser.search("-nMC"))
        nMC = parser.next(nMC);

    // ODE
    odeDef ODE;
    ODE.ode = FITZNAG;
    setProblem(&ODE);

    // Reference solution
    Butcher refTableau(RK4, EXPLICIT, 0);
    RungeKutta refSolver(ODE, refTableau);
    double hRef = 1e-5;
    auto NRef = static_cast<unsigned int>(std::round(T / hRef));
    VectorXd refSolution = ODE.initialCond;
    double refPhi;
    for (unsigned int j = 0; j < NRef; j++) {
        refSolution = refSolver.oneStep(hRef, refSolution, ODE.refParam);
    }
    refPhi = phi(refSolution);

    // Error computation
    Butcher tableau(EXPTRAPEZ, EXPLICIT, 0);
    std::random_device device;
    std::default_random_engine generator{device()};

    std::ofstream output(DATA_PATH + outputFileName + ".txt", std::ofstream::out | std::ofstream::trunc);
    for (int k = 0; k < nExp; k++) {
        std::cout << "Computation for h = " << h << std::endl;
        auto N = static_cast<unsigned int>(std::round(T / h));
        double expPhi = 0.0;

        int j = 0;
        #pragma omp parallel for num_threads(6) reduction(+:expPhi) private(j)
        for (j = 0; j < nMC; j++) {
            RungeKuttaRandomH probSolver(&generator, ODE, tableau, h, p);
            VectorXd solution = ODE.initialCond;
            for (unsigned int i = 0; i < N; i++) {
                solution = probSolver.oneStep(solution, ODE.refParam);
            }
            expPhi += phi(solution);
        }
        double error = std::abs(expPhi / nMC - refPhi);

        output << h << "\t" << error << std::endl;
        h /= 2.0;
    }
    output.close();

    return 0;
}
