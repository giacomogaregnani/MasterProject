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
    int nExp = 1, nMC = 10, nMCExt = 10;

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
    if (parser.search("-nMCExt"))
        nMCExt = parser.next(nMCExt);

    // ODE
    odeDef ODE;
    ODE.ode = FITZNAG;
    setProblem(&ODE);

    // Reference solution
    Butcher refTableau(RK4, EXPLICIT, 0);
    RungeKutta refSolver(ODE, refTableau);
    double hRef = 1e-6;
    auto NRef = static_cast<unsigned int>(std::round(T / hRef));
    VectorXd refSolution = ODE.initialCond;

    for (unsigned int j = 0; j < NRef; j++) {
        refSolution = refSolver.oneStep(hRef, refSolution, ODE.refParam);
    }
    double refPhi = phi(refSolution);

    // Error computation
    Butcher tableau(RK4, EXPLICIT);
    std::random_device device;
    std::default_random_engine generator{device()};
    auto N = static_cast<unsigned int>(std::round(T / h));

    std::ofstream output(DATA_PATH + outputFileName + ".txt", std::ofstream::out | std::ofstream::trunc);
    for (int k = 0; k < nExp; k++) {

        std::cout << "Computation for M = " << nMC << std::endl;
        double error = 0.0;

        for (int l = 0; l < nMCExt; l++) {
            double est = 0.0;
            int j = 0;
            for (j = 0; j < nMC; j++) {
                RungeKuttaRandomH probSolver(&generator, ODE, tableau, h, p);
                VectorXd solution = ODE.initialCond;
                for (unsigned int i = 0; i < N; i++) {
                    solution = probSolver.oneStep(solution, ODE.refParam);
                }
                est += phi(solution);
            }
            est /= nMC;
            error += (est - refPhi) * (est - refPhi);
        }
        error = error / nMCExt;

        output << nMC << "\t" << sqrt(error) << std::endl;
        nMC *= 2;
    }
    output.close();

    return 0;
}
