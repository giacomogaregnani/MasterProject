#include <RandomTimeStep.hpp>
#include <RandomCalibration.hpp>
#include <fstream>
#include <iomanip>
#include <GetPot>

void mAndS(std::vector<VectorXd>& data, VectorXd& m, MatrixXd& s)
{
    m = VectorXd::Zero(m.size());
    for (const auto &it : data)
        m += it;

    m /= data.size();

    s = MatrixXd::Zero(m.size(), m.size());

    for (const auto &it : data)
        s += (it - m) * (it - m).transpose();
    s /= (data.size() - 1);
}

double hamKepler(VectorXd& v, double h = 0)
{
    return 0.5 * (v(0) * v(0) + v(1) * v(1)) - 1.0 / (sqrt(v(2) * v(2) + v(3) * v(3)));
}

double hamPendulum(VectorXd& v, double h = 0)
{
    return v(0) * v(0) / 2.0 - cos(v(1));
}

int main(int argc, char* argv[])
{
    GetPot parser(argc, argv);

    double h = 0.05, p = 1.5,
            T = 10, sigma = 0.01;
    int nMC = 10, nMCMC = 100;
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
    if (parser.search("-nMCMC"))
        nMCMC = parser.next(nMCMC);
    if (parser.search("-sigma"))
        sigma = parser.next(sigma);
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
            case 4:
                whichODE = KEPLER;
                break;
            case 5:
                whichODE = PENDULUM;
                break;
            case 6:
                whichODE = HIRES;
                break;
            case 7:
                whichODE = TEST1D;
                break;
            case 8:
                whichODE = MODIFPEND;
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
    if (whichODE == TEST1D) {
        ODE.refParam[0] = -0.1;
    }

    // Integration parameters
    auto N = static_cast<unsigned int>(T / h);

    // Random seed
    std::default_random_engine generator{(unsigned int) time(nullptr)};

    // Output files
    std::string outputName = std::string(DATA_PATH) + outputFile + ".txt";
    std::string outputName2 = std::string(DATA_PATH) + outputFile + "Growth.txt";
    std::string outputName3 = std::string(DATA_PATH) + outputFile + "Growth2.txt";
    std::ofstream output(outputName, std::ofstream::out | std::ofstream::trunc);
    std::ofstream output2(outputName2, std::ofstream::out | std::ofstream::trunc);
    std::ofstream output3(outputName3, std::ofstream::out | std::ofstream::trunc);

    // Numerical method
    Butcher tableau(EXPTRAPEZ, EXPLICIT);
    Butcher embeddedTableau(GAUSS6, IMPLICIT);
    Butcher highOrderTableau(GAUSS6, IMPLICIT);

    RungeKutta detMethod(ODE, tableau);
    RungeKutta detMethod2(ODE, embeddedTableau);
    RungeKutta refMethod(ODE, highOrderTableau);
    RungeKuttaRandomH probMethod(&generator, ODE, tableau, h, p);

    // Calibration phase
    RandomCalibration<RungeKuttaRandomH> calibrator(h*200, &probMethod, &detMethod, &detMethod2);
    calibrator.calibrate(nMC, nMCMC, sigma);

    std::vector<double> constants = calibrator.getConstants();
    std::vector<double> densities = calibrator.getDensities();

    // Choose the MAP
    auto maxDensity = std::max_element(densities.begin(), densities.end());
    auto indexMax = std::distance(densities.begin(), maxDensity);
    double calibrationConstant = constants[indexMax];
    probMethod.setConstant(calibrationConstant);
    std::cout << calibrationConstant << std::endl;

    for (int i = 0; i < constants.size(); i++) {
        output << constants[i] << "\t" << densities[i] << std::endl;
    }

    // Reference solution
    std::cout << "Generating the reference solution..." << std::endl;
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

    // Deterministic integration (estimate the error)
    std::cout << "Error estimation..." << std::endl;
    VectorXd detSolution = ODE.initialCond;
    double error;
    for (unsigned int i = 0; i < N; i++) {
        error = (detSolution - refSolution[i]).norm();
        detSolution = detMethod.oneStep(h, detSolution, ODE.refParam);
        output2 << (h*i) << "\t" << error << std::endl;
    }

    // Probabilistic calibrated integration (with square root of time)
    std::cout << "Probabilistic integrator, square root of time..." << std::endl;
    std::vector<VectorXd> probSolutions(nMC, ODE.initialCond);
    MatrixXd var(ODE.size, ODE.size); VectorXd mean(ODE.size);
    double C = calibrationConstant;
    for (unsigned int i = 0; i < N; i++) {
        if (h * i > 1) {
            C = calibrationConstant * std::pow(h * i, 0.5);
        }
        probMethod.setConstant(C);
        mAndS(probSolutions, mean, var);
        output2 << (h*i) << "\t" << std::sqrt(var.norm()) << std::endl;
        for (unsigned int j = 0; j < nMC; j++) {
            probSolutions[j] = probMethod.oneStep(probSolutions[j], ODE.refParam);
        }
    }

    // Probabilistic calibrated integration (without square root of time)
    std::cout << "Probabilistic integrator, fixed..." << std::endl;
    probMethod.setConstant(calibrationConstant);
    for (unsigned int i = 0; i < nMC; i++)
        probSolutions[i] = ODE.initialCond;

    double mHamiltonian = hamPendulum(ODE.initialCond);

    for (unsigned int i = 0; i < N; i++) {
        mAndS(probSolutions, mean, var);
        output3 << (h*i) << "\t" << std::sqrt(var.norm()) << "\t" << mHamiltonian << std::endl;

        mHamiltonian = 0.0;
        for (unsigned int j = 0; j < nMC; j++) {
            probSolutions[j] = probMethod.oneStep(probSolutions[j], ODE.refParam);
            mHamiltonian += hamPendulum(probSolutions[j]);
        }
        mHamiltonian /= nMC;
    }

    output.close();
    output2.close();
    output3.close();

    return 0;
}