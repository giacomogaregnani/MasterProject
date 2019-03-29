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
    if (parser.search("-nMC"))
        nMC = parser.next(nMC);
    if (parser.search("-tol"))
        tol = parser.next(tol);
    if (parser.search("-nExp"))
        nExp = parser.next(nExp);
    if (parser.search("-p"))
        p = parser.next(p);

    // ODE
    odeDef ODE;
    ODE.ode = VDPOL;
    setProblem(&ODE);
    ODE.refParam[0] = 1;

    // Computation
    std::default_random_engine generator{(unsigned int) time(nullptr)};

    Butcher tableau(EXPTRAPEZ, EXPLICIT);
    Butcher embTableau(EULERFORWARD, EXPLICIT);
    RungeKutta embSolverOne(ODE, tableau);
    RungeKutta embSolverTwo(ODE, embTableau);

    // Reference solution
    Butcher refTableau(RK4, EXPLICIT);
    RungeKutta refSolver(ODE, refTableau);
    double refH = 1e-6;
    auto refN = static_cast<unsigned int>(T / refH);

    VectorXd refSolution = ODE.initialCond;
    for (unsigned int i = 0; i < refN; i++) {
        refSolution = refSolver.oneStep(refH, refSolution, ODE.refParam);
    }

    // Adaptive RTS-RK
    std::ofstream summary(DATA_PATH + outputFileName + "summary.txt", std::ofstream::out | std::ofstream::trunc);

    for (unsigned int k = 0; k < nExp; k++) {

        std::cout << std::fixed << std::setprecision(5) << "tolerance = " << tol << std::endl;

        // EMBEDDED METHODS
        VectorXd embSolutionOne = ODE.initialCond;
        VectorXd embSolutionTwo = ODE.initialCond;
        VectorXd tmpSolution;
        double errEst = 0.0, oldErrEst = 0.0;
        int nTotEmb = 0, nRejEmb = 0;
        double meanHEmb = 0.0, hEmb = h, hOldEmb = h, hNewEmb, fac, fac1, fac2, fac3, tolEmb = std::sqrt(tol);

        std::ofstream outputEmb(DATA_PATH + outputFileName + std::to_string(k) + "Emb.txt", std::ofstream::out | std::ofstream::trunc);
        double tEmb = 0;

        double locConst = 0.0;

        while (tEmb < T) {
            std::cout << std::fixed << std::setprecision(5) << "t = " << tEmb << "\t"
                      << "h = " << hEmb << "\t";

            if (tEmb + hEmb > T)
                hEmb = T - tEmb;

            tmpSolution = embSolutionOne;
            embSolutionOne = embSolverOne.oneStep(hEmb, embSolutionOne, ODE.refParam);
            embSolutionTwo = embSolverTwo.oneStep(hEmb, embSolutionOne, ODE.refParam);

            errEst = (embSolutionOne - embSolutionTwo).norm();
            locConst += errEst / hEmb;

            if (tEmb == 0) {
                oldErrEst = errEst;
                hOldEmb = hEmb;
            }

            std::cout << "err = " << errEst << "\t";

            fac1 = std::sqrt(tolEmb / errEst);
            fac2 = std::sqrt(oldErrEst / errEst);
            fac3 = hEmb / hOldEmb;
            fac = fac1;//  * fac2 * fac3;
            hNewEmb = std::min(5.0, std::max(0.1,  0.9 * fac)) * hEmb;
            std::cout << "fac = " << fac1 << " " << fac2 << " " << fac3 << " " << fac << "\t";

            if (errEst <= tolEmb) {
                tEmb += hEmb;
                meanHEmb += hEmb;
                outputEmb << tEmb << "\t" << embSolutionOne.transpose() << "\t" << errEst << std::endl;
            } else {
                nRejEmb++;
                embSolutionOne = tmpSolution;
                std::cout << "rejected" << "\t";
            }

            oldErrEst = errEst;
            hOldEmb = hEmb;
            hEmb = hNewEmb;

            nTotEmb++;

            std::cout << std::endl;
        }

        double embErr = (embSolutionOne - refSolution).norm();
        summary << tol << "\t" << embErr << "\t" << nTotEmb * tableau.getStages() << "\t"
                << static_cast<double>(nRejEmb) / nTotEmb << "\t" << meanHEmb / (nTotEmb - nRejEmb) << std::endl;

        locConst /= nTotEmb;
        RungeKuttaAddNoise probSolver(&generator, ODE, tableau, h, p, locConst);

        // PROBABILISTIC ADAPTATION
        std::vector<VectorXd> solution(nMC, ODE.initialCond);
        std::vector<VectorXd> newSolution = solution;
        VectorXd mean = VectorXd::Zero(ODE.size);

        // Output file
        std::ofstream output(DATA_PATH + outputFileName + std::to_string(k) + ".txt", std::ofstream::out | std::ofstream::trunc);
        double t = 0, stddev = 0, oldStd = 0, controller = 0, hProb = h, hOldProb = h, hNewProb, hMean = 0;
        int nRej = 0, nTot = 0;

        probSolver.setH(h);

        while (t < T) {

            std::cout << std::fixed << std::setprecision(15) << "t = " << t << "\t"
                      << "h = " << hProb << "\t";

            if (t + hProb > T) {
                hProb = T - t;
                probSolver.setH(hProb);
            }

            for (int i = 0; i < nMC; i++)
                newSolution[i] = probSolver.oneStep(solution[i], ODE.refParam);

            stddev = computeStd(newSolution);

            std::cout << "oldStd = " << oldStd
                      << "\t newStd = " << stddev
                      << "\t newStd - oldStd = " << stddev - oldStd << "\t";

            if (stddev - oldStd > 0) {
                errEst = std::pow(stddev - oldStd, 1.0 / 1.5);
            } else {
                errEst = 0.0;
            }

            if (t == 0) {
                oldErrEst = errEst;
                hOldProb = hProb;
            }

            fac1 = std::sqrt(tol / errEst);
            fac2 = std::sqrt(oldErrEst / errEst);
            fac3 = hProb / hOldProb;
            controller = fac1; // * fac2 * fac3;
            hNewProb = std::min(5.0, std::max(0.1, 0.9 * controller)) * hProb;

            if (errEst <= tol) {
                t += hProb;
                output << t << "\t";
                for (auto it : solution)
                    output << it.transpose() << "\t";
                output << errEst << std::endl;
                solution = newSolution;
                oldStd = stddev;
                hMean += hProb;
            } else {
                nRej++;
                std::cout << "rejected" << "\t";
            }

            nTot++;

            hOldProb = hProb;
            oldErrEst = errEst;
            hProb = hNewProb;
            probSolver.setH(hProb);

            std::cout << std::endl;

        }
        output.close();

        double err = 0;
        VectorXd meanSol = VectorXd::Zero(refSolution.size());
        for (const auto &it : newSolution) {
            meanSol += it;
            err += (it - refSolution).dot(it - refSolution);
        }
        meanSol /= nMC;
        double meanErr = (meanSol - refSolution).norm();
        err = std::sqrt(err / nMC);

        summary << tol << "\t" << err << "\t" << nTot * tableau.getStages() << "\t"
                << static_cast<double>(nRej) / nTot << "\t" << hMean / (nTot - nRej) << "\t" << meanErr << "\t" << stddev << std::endl;


        std::cout << "rejection ratio = " << static_cast<double>(nRej) / nTot << std::endl;
        std::cout << "loc const = " << locConst << std::endl;
        tol /= 2;
    }


    return 0;
}