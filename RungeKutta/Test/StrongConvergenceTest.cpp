#include <RandomTimeStep.hpp>
#include <UtilitiesGG.hpp>
#include <iomanip>

int main(int argc, char* argv[])
{
    odeDef ODE;
    ODE.ode = FITZNAG;
    setProblem(&ODE);

    std::vector<double> param = {0.2, 0.2, 3.0};

    Butcher tableau(RK4, EXPLICIT, 0);

    double T = 1.0, h = 0.01;
    unsigned int nMC = 10000, nExp = 5;

    std::default_random_engine generator{(unsigned int) time(NULL)};

    // Output file
    std::string outputFilename = std::string(DATA_PATH) + argv[1] + ".txt";
    std::ofstream output;
    output.open(outputFilename, std::ofstream::out | std::ofstream::trunc);

    // Reference solution
    VectorXd refSolution(ODE.size);
    if (ODE.ode == TEST1D) {
        refSolution(0) = exp(param[0] * T);
    } else {
        std::string refFilename = std::string(DATA_PATH) + argv[2] + ".txt";
        refSolution = loadRefSolution(refFilename, ODE.size);
    }

    std::vector<double> errors(nExp);
    std::vector<double> orders(nExp - 1);

    // "Timestep" loop
    for (unsigned int j = 0; j < nExp; j++) {

        VectorXd MCMean = VectorXd::Zero(ODE.size);
        unsigned int N = static_cast<unsigned int>(T / h);

        std::cout << "Mean timestep = "
                  << h
                  << "\t\t"
                  << "Number of steps = "
                  << N
                  << std::endl;

        // MC loop
        unsigned int k;
        std::vector<double> finalTime(nMC), MCError(nMC);

        #pragma omp parallel for num_threads(30) private(k)
        for (k = 0; k < nMC; k++) {

            RungeKuttaRandomH Method(&generator, ODE, param, tableau, h, std::atof(argv[3]));
            VectorXd solution = ODE.initialCond;

            // Time integration loop
            for (unsigned int i = 0; i < N; i++) {
                solution = Method.oneStep(solution);
            }

            MCError[k] = (solution - refSolution).norm();
            finalTime[k] = Method.getCurrentTime();
        }

        double strongError = 0.0;
        for (auto it : MCError) {
            strongError += it;
        }
        strongError /= nMC;

        double meanFinalTime = 0.0;
        for (auto it : finalTime) {
            meanFinalTime += it;
        }
        meanFinalTime /= nMC;
        std::cout << "Mean final time = " << meanFinalTime << std::endl;

        errors[j] = strongError;
        output << std::fixed << std::setprecision(20) << h << "\t" << strongError << "\n";

        h /= 2;
    }

    computeOrderOfConvergence(errors, orders, 2);
    printConvergenceInfo(errors, orders);

    output.close();

    return 0;
}
