#include <RandomTimeStep.hpp>
#include <UtilitiesGG.hpp>
#include <iomanip>

int main(int argc, char* argv[])
{
    if (argc < 6) {
        throw std::invalid_argument("Inputs:\n"
                                    "h\n"
                                    "T\n"
                                    "nMC\n"
                                    "p\n"
                                    "output file\n");
    }

    odeDef ODE;
    ODE.ode = LORENZ;
    setProblem(&ODE);

    std::vector<double> param = ODE.refParam;

    Butcher tableau(RK4, EXPLICIT, 0);

    // Parse input
    double T = std::atof(argv[2]), h = std::atof(argv[1]);
    unsigned long nMC = std::stoul(argv[3]);
    double p = std::atof(argv[4]);
    unsigned int N = static_cast<unsigned int>(T / h);

    std::default_random_engine generator{(unsigned int) time(NULL)};

    // Output file
    std::string outputFilenameInfo = std::string(DATA_PATH) + "tmpTraj/" + argv[5] + "Info.txt";
    std::string outputFilename = std::string(DATA_PATH) + "tmpTraj/" + argv[5] + ".txt";
    std::ofstream output, infoFile;
    output.open(outputFilename, std::ofstream::out | std::ofstream::trunc);

    // Write on Info File
    infoFile.open(outputFilenameInfo, std::ofstream::out | std::ofstream::trunc);
    infoFile << h << "\t" << T << "\t" << nMC;
    infoFile.close();


    // Initialize
    RungeKuttaRandomH Method(&generator, ODE, param, tableau, h, p);

    for (unsigned long k = 0; k < nMC; k++) {

        VectorXd solution = ODE.initialCond;
        output << solution.transpose() << std::endl;

        // Time integration loop
        for (unsigned int i = 0; i < N; i++) {
            solution = Method.oneStep(solution);
            output << solution.transpose() << std::endl;
        }

    }

    output.close();

    return 0;
}
