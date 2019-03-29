#include <RandomTimeStep.hpp>
#include <fstream>
#include <iomanip>
#include <GetPot>
#include <queue>

double computeStd(VectorXd& mean, std::vector<VectorXd> &data)
{
    mean = VectorXd::Zero(mean.size());
    for (const auto &it : data)
        mean += it;

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
    int nMC = 10;
    problems whichODE = FITZNAG;

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
    if (parser.search("-problem")) {
        switch (parser.next(0)){
            case 0:
                whichODE = FITZNAG;
                break;
            case 1:
                whichODE = HENHEIL;
                break;
            case 2:
                whichODE = LORENZ;
                break;
            case 3:
                whichODE = VDPOL;
                break;
        }
    }

    // ODE
    odeDef ODE;
    ODE.ode = whichODE;
    setProblem(&ODE);
    if (whichODE == VDPOL) {
        ODE.refParam[0] = 1;
    }

    // Numerical method
    Butcher tableau(EXPTRAPEZ, EXPLICIT);
    Butcher highOrderTableau(RK4, EXPLICIT);

    // Integration parameters
    auto N = static_cast<unsigned int>(T / h);

    // Random seed
    std::default_random_engine generator{(unsigned int) time(NULL)};

    // Output files
    std::string outputProbName = std::string(DATA_PATH) + outputFile + ".txt";
    std::string outputDetName = std::string(DATA_PATH) + outputFile + "det.txt";
    std::ofstream outputDet(outputDetName, std::ofstream::out | std::ofstream::trunc);
    std::ofstream outputProb(outputProbName, std::ofstream::out | std::ofstream::trunc);

    // Reference solution
    RungeKutta refMethod(ODE, highOrderTableau);
    double hRef = h / 100;
    unsigned int NRef = N * 100;
    std::vector<VectorXd> refSolution(N, ODE.initialCond);
    VectorXd tempSolution = refSolution[0];
    int k = 0;
    for (unsigned int j = 0; j  < NRef; j++) {
        if (j % 100 == 0) {
            refSolution[k++] = tempSolution;
        }
        tempSolution = refMethod.oneStep(hRef, tempSolution, ODE.refParam);
    }

    // Deterministic
    RungeKutta detMethod(ODE, tableau);
    VectorXd solution = ODE.initialCond;
    VectorXd solution2 = solution;
    VectorXd temp = ODE.initialCond;
    double errConstant = 0;

    for (unsigned int j = 0; j < N; j++) {
        outputDet << std::fixed << std::setprecision(10)
                  << h * j << "\t"
                  << (solution - refSolution[j]).norm() << std::endl;

        temp = solution;
        solution = detMethod.oneStep(h, solution, ODE.refParam);
        solution2 = refMethod.oneStep(h, temp, ODE.refParam);

        // errConstant += (solution - solution2).norm();
        errConstant = std::max(errConstant, (solution - solution2).norm());
    }
    errConstant /= std::pow(h, p+1);

    // double noiseConstant = errConstant / (std::pow(h, p+1) * N);
    double noiseConstant = errConstant * std::sqrt(h);
    std::cout << noiseConstant << std::endl;
    outputDet.close();

    // Initialization
    RungeKuttaAddNoise probMethod(&generator, ODE, tableau, h, p, noiseConstant);
    std::vector<VectorXd> probSolutions(nMC, ODE.initialCond);
    double s = 0.0;
    VectorXd mean(ODE.initialCond);

    for (unsigned int j = 0; j < N; j++) {
        outputProb << std::fixed << std::setprecision(10)
                   << h * j << "\t"
                   << s << std::endl;

        for (unsigned int i = 0; i < nMC; i++)
            probSolutions[i] = probMethod.oneStep(probSolutions[i], ODE.refParam);

        s = computeStd(mean, probSolutions);
    }

    outputProb.close();

    return 0;
}