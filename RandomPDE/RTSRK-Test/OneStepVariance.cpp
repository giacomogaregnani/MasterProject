#include <RandomTimeStep.hpp>
#include <fstream>
#include <iomanip>
#include <GetPot>

double computeStd(std::vector<VectorXd> &data)
{
    VectorXd mean = VectorXd::Zero(data[0].size());
    for (const auto &it : data) {
        mean += it;
    }
    mean /= data.size();

    double varNorm = 0.0;

    for (const auto &it : data)
        varNorm += (it - mean).dot(it - mean);
    varNorm /= (data.size() - 1);

    return std::sqrt(varNorm);
}

int main(int argc, char* argv[])
{
    GetPot parser(argc, argv);

    double h = 0.1, p = 1.5,
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
    if (parser.search("-p"))
        p = parser.next(p);
    if (parser.search("-nExp"))
        nExp = parser.next(nExp);

    // ODE
    odeDef ODE;
    ODE.ode = FITZNAG;
    setProblem(&ODE);

    // Numerical method
    Butcher tableau(EULERFORWARD, EXPLICIT);

    // Integration parameters
    auto N = static_cast<unsigned int>(T / h);

    // Random seed
    std::default_random_engine generator{(unsigned int) time(NULL)};
    RungeKuttaAddNoise Method(&generator, ODE, tableau, h, p);

    for (int k = 0; k < nExp; k++) {

        std::string outputFileName = std::string(DATA_PATH) + outputFile + std::to_string(k) + ".txt";
        std::ofstream output(outputFileName, std::ofstream::out | std::ofstream::trunc);
        std::vector<VectorXd> solution(nMC, ODE.initialCond);

        double std = 0.0;

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < nMC; j++)
                solution[j] = Method.oneStep(solution[j], ODE.refParam);

            std = computeStd(solution);
            output << h * (i + 1) << "\t" << std << std::endl;
        }

        output.close();

        h /= 2;
        N *= 2;
        Method.setH(h);

    }

    return 0;
}