#include <RandomTimeStep.hpp>
#include <fstream>
#include <iomanip>
#include <GetPot>
#include <iostream>

int main(int argc, char* argv[])
{
    GetPot parser(argc, argv);

    double h = 0.1,
            T = 1,
            p = 1.0;
    int nMC = 10;

    std::string outputFileName;

    if (parser.search("-h"))
        h = parser.next(h);
    if (parser.search("-T"))
        T = parser.next(T);
    if (parser.search("-output"))
        outputFileName = parser.next(" ");
    if (parser.search("-p"))
        p = parser.next(p);
    if (parser.search("-nMC"))
        nMC = parser.next(nMC);

    // ODE
    odeDef ODE;
    ODE.ode = PEROX;
    setProblem(&ODE);

    // Error computation
    Butcher tableau(EXPTRAPEZ, EXPLICIT, 0);
    // Butcher tableau(RKC, STABEXP, 10);
    std::default_random_engine generator{(unsigned int) time(NULL)};

    std::ofstream output(DATA_PATH + outputFileName + ".txt", std::ofstream::out | std::ofstream::trunc);

    std::cout << "Computation for h = " << h << std::endl;
    unsigned int N = static_cast<unsigned int>(std::round(T / h));

    RungeKuttaRandomH probSolver(&generator, ODE, tableau, h, p);
    RungeKutta detSolver(ODE, tableau);

    std::vector<VectorXd> solution(nMC, ODE.initialCond);
    VectorXd detSolution = ODE.initialCond;

    for (unsigned int i = 0; i < N; i++) {
        detSolution = detSolver.oneStep(h, detSolution, ODE.refParam);

        for (int j = 0; j < nMC; j++) {
            solution[j] = probSolver.oneStep(solution[j], ODE.refParam);
        }

        if (i % 50 == 0) {
            double error = 0.0;
            for (int j = 0; j < nMC; j++) {
                error += sqrt((detSolution - solution[j]).dot(detSolution - solution[j]));
            }
            output << h*i << "\t" << error / nMC << std::endl;
        }
    }

    output.close();

    return 0;
}