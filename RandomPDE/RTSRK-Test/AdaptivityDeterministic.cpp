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

    MatrixXd var = MatrixXd::Zero(mean.size(), mean.size());
    for (const auto &it : data) {
        var += (it - mean)*(it - mean).transpose();
    }
    var /= (data.size() - 1);

    return std::sqrt(var.norm());
}

int main(int argc, char* argv[])
{
    GetPot parser(argc, argv);

    double h = 0.1,
            T = 1,
            tol = 1.0,
            p = 2.5;
    int nMC = 10, nExp = 1;

    std::string outputFileName;

    if (parser.search("-h"))
        h = parser.next(h);
    if (parser.search("-T"))
        T = parser.next(T);
    if (parser.search("-output"))
        outputFileName = parser.next(" ");
    if (parser.search("-tol"))
        tol = parser.next(tol);
    if (parser.search("-nExp"))
        nExp = parser.next(nExp);

    // ODE
    odeDef ODE;
    ODE.ode = LORENZ;
    setProblem(&ODE);

    // Computation
    Butcher tableau(EXPTRAPEZ, EXPLICIT);
    std::default_random_engine generator{(unsigned int) time(nullptr)};
    RungeKutta embSolverOne(ODE, tableau);
    Butcher embTableau(EULERFORWARD, EXPLICIT);
    RungeKutta embSolverTwo(ODE, embTableau);

    // Reference solution
    /* Butcher refTableau(RK4, EXPLICIT);
    RungeKutta refSolver(ODE, refTableau);
    double refH = 1e-4;
    auto refN = static_cast<unsigned int>(T / refH);

    VectorXd refSolution = ODE.initialCond;
    for (unsigned int i = 0; i < refN; i++) {
        refSolution = refSolver.oneStep(refH, refSolution, ODE.refParam);
    } */

    std::vector<double> hZeros(nExp);
    std::uniform_real_distribution<double> hZeroDistribution(h - h * h, h + h * h);

    for (unsigned int k = 0; k < nExp; k++) {
        std::cout << std::fixed << std::setprecision(5) << "tolerance = " << tol << std::endl;

        // EMBEDDED METHODS (DO NOT TOUCH THE CODE WORKS PERFECTLY FINE)
        VectorXd embSolutionOne = ODE.initialCond;
        VectorXd embSolutionTwo = ODE.initialCond;
        VectorXd tmpSolution;
        double errEst = 0.0, oldErrEst = 0.0;
        int nTotEmb = 0, nRejEmb = 0;
        double hEmb = hZeroDistribution(generator);
        double hOldEmb = h, hNewEmb;

        std::ofstream outputEmb(DATA_PATH + outputFileName + std::to_string(k) + ".txt",
                                std::ofstream::out | std::ofstream::trunc);
        double tEmb = 0;

        while (tEmb < T) {
            std::cout << std::fixed << std::setprecision(5) << "t = " << tEmb << "\t"
                      << "h = " << hEmb << "\t";

            if (tEmb + hEmb > T)
                hEmb = T - tEmb;

            tmpSolution = embSolutionOne;
            embSolutionOne = embSolverOne.oneStep(hEmb, embSolutionOne, ODE.refParam);
            embSolutionTwo = embSolverTwo.oneStep(hEmb, embSolutionOne, ODE.refParam);

            VectorXd d = embSolutionOne - embSolutionTwo;
            errEst = std::sqrt(d.dot(d));
            if (tEmb == 0)
                oldErrEst = errEst;

            std::cout << "err = " << errEst << "\t";

            nTotEmb++;
            double fac1 = std::sqrt(tol / errEst);
            double fac2 = std::sqrt(oldErrEst / errEst);
            double fac3 = hEmb / hOldEmb;
            double fac = fac1 * fac2 * fac3;
            hNewEmb = 0.9 * hEmb * std::min(5.0, std::max(0.1, fac));
            std::cout << "fac = " << fac1 << " " << fac2 << " " << fac3 << " " << fac << "\t";

            if (errEst <= tol) {
                tEmb += hEmb;
                outputEmb << tEmb << "\t" << embSolutionOne.transpose() << std::endl;
            } else {
                nRejEmb++;
                embSolutionOne = tmpSolution;
                std::cout << "rejected" << "\t";
            }

            oldErrEst = errEst;
            hOldEmb = hEmb;
            hEmb = hNewEmb;

            std::cout << std::endl;
        }
    }

    return 0;
}
