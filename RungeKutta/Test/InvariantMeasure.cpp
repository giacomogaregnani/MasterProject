#include <RandomTimeStep.hpp>
#include <UtilitiesGG.hpp>
#include <iomanip>
#include <GetPot>

int main(int argc, char* argv[])
{
    GetPot parser(argc, argv);

    double h = 0.01;
    std::string outputName = "output";
    int p = 1;
    double T = 50;
    int nM = 1;

    if (parser.search("-h"))
        h = parser.next(h);
    if (parser.search("-output"))
        outputName = parser.next("output");
    if (parser.search("-p"))
        p = parser.next(p);
    if (parser.search("-T"))
        T = parser.next(T);
    if (parser.search("-nM"))
        nM = parser.next(nM);

    odeDef ODE;
    ODE.ode = LORENZ;
    setProblem(&ODE);
    int N = static_cast<int> (std::round(T / h));

    Butcher tableau(IMPMID, IMPLICIT, 0);

    std::default_random_engine generator{(unsigned int) time(NULL)};

    // Output file
    std::string outputFilename = std::string(DATA_PATH) + "Measure/" + outputName + ".txt";
    std::ofstream output;
    output.open(outputFilename, std::ofstream::out | std::ofstream::trunc);


    // Time average
    /* RungeKutta detSolver(ODE, ODE.refParam, tableau);
    VectorXd detSolution = ODE.initialCond;
    double detMean = 0;
    double hDet = 1e-4;
    int detN = static_cast<int> (std::round(200.0 / hDet));

    for (int i = 0; i < detN; i++) {
        detSolution = detSolver.oneStep(detSolution, hDet);
        detMean += detSolution(1);
    }
    detMean /= detN;
    output << 0 << "\t" << std::fixed << std::setprecision(20) << detMean << std::endl;
    */

    int maxM = static_cast<int>(std::round(pow(2.0, nM)));
    // Space averages
    for (int M = 1; M < maxM + 1; M *= 2) {

        double mean = 0;
        int k;
        std::vector<double> results(static_cast<unsigned long>(M));

        #pragma omp parallel for num_threads(20) private(k)
        for (k = 0; k < M; k++) {
            RungeKuttaRandomH probSolver(&generator, ODE, ODE.refParam, tableau, h, p);
            VectorXd solution = ODE.initialCond;

            // Time integration loop
            for (int i = 0; i < N; i++) {
                solution = probSolver.oneStep(solution);
            }

            results[k] = solution(1);
        }

        for (auto it : results) {
            mean += it;
        }

        mean /= M;
        output << M << "\t" << std::fixed << std::setprecision(20) << mean << std::endl;

        if (M == maxM) {
            std::string fullFileName = std::string(DATA_PATH) + "Measure/" + outputName + "_full.txt";
            std::fstream fullOutput(fullFileName, std::ofstream::out | std::ofstream::trunc);
            for (auto it : results) {
                fullOutput << it << std::endl;
            }
            fullOutput.close();
        }
    }

    output.close();

    return 0;
}
