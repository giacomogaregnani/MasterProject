#include "RandomCalibration.hpp"

double bDistanceScalar(double m1, double v1, double m2, double v2)
{
    return 0.25 * std::log(0.25 * (v1 / v2 + v2 / v1 + 2.0)) + 0.25 * ((m1 - m2) * (m1 - m2) / (v1 + v2));
}

void computeMeanAndVariance(std::vector<VectorXd>& data, VectorXd& m, MatrixXd& s)
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

void computeMeanAndVariance(std::vector<double>& data, double& m, double& s)
{
    m = 0.0;
    for (const auto &it : data)
        m += it;

    m /= data.size();

    s = 0.0;

    for (const auto &it : data)
        s += (it - m) * (it - m);
    s /= (data.size() - 1);
}

template <class T>
RandomCalibration<T>::RandomCalibration(double finalTime, T* method,
                                        RungeKutta* detMethod, RungeKutta* embMethod):
        finalTime(finalTime),
        probMethod(method),
        detMethod(detMethod),
        embMethod(embMethod)
{
    h = probMethod->getH();
    ODE = probMethod->getODE();
}

template <class T>
void RandomCalibration<T>::embErrorsEstimate()
{
    auto N = static_cast<unsigned int>(std::round(finalTime / h));
    VectorXd solution = ODE.initialCond;
    VectorXd embSolution = ODE.initialCond;
    VectorXd temp(solution.size());

    for (unsigned int i = 0; i < N; i++) {
        temp = solution;
        solution = detMethod->oneStep(h, solution, ODE.refParam);
        embSolution = embMethod->oneStep(h, embSolution, ODE.refParam);
        detSolution.push_back(solution);
        embErrors.push_back(solution - embSolution);
    }
}

template <class T>
double RandomCalibration<T>::distribution(double C, unsigned int nMC)
{
    probMethod->setConstant(C);

    auto N = static_cast<unsigned int>(std::round(finalTime / h));
    std::vector<VectorXd> solutions(nMC, ODE.initialCond);
    double density = 0.0;
    VectorXd solutionMean(ODE.size);
    MatrixXd solutionVariance(ODE.size, ODE.size);
    MatrixXd embErrorMat;

    for (unsigned int i = 0; i < N; i++) {
        for (unsigned int j = 0; j < nMC; j++) {
            solutions[j] = probMethod->oneStep(solutions[j], ODE.refParam);
        }
        computeMeanAndVariance(solutions, solutionMean, solutionVariance);

        double distance = 1;
        for (unsigned int k = 0; k < ODE.size; k++) {
            distance *= bDistanceScalar(detSolution[i](k), embErrors[i](k) * embErrors[i](k),
                                        solutionMean(k), solutionVariance(k, k));
        }
        density -= distance;
    }

    return density;
}

template <class T>
void RandomCalibration<T>::calibrate(unsigned int nMC, unsigned int nMCMC, double sigma)
{
    embErrorsEstimate();

    double alpha;
    unsigned int accRatio = 0;

    samplesConstant.push_back(1.0);
    densitiesConstant.push_back(distribution(samplesConstant.back(), nMC));
    double proposedValue, proposedDensity;

    std::default_random_engine generator{(unsigned int) time(nullptr)};
    std::normal_distribution<double> proposal(0, sigma);
    std::uniform_real_distribution<double> probGen(0.0, 1.0);

    for (int i = 0; i < nMCMC; i++) {
        if (i % 50 == 0) {
            std::cout << "Completed " << i << " iterations out of " << nMCMC << std::endl
                      << "Current acc. ratio = " << 1.0 / (i + 1) * accRatio << std::endl;
            std::cout << samplesConstant.back() << "\t" << densitiesConstant.back() << std::endl;
        }

        proposedValue = samplesConstant.back() + proposal(generator);
        proposedDensity = distribution(proposedValue, nMC);

        alpha = proposedDensity - densitiesConstant.back();
        alpha = std::min(0.0, alpha);

        if (std::log(probGen(generator)) < alpha) {
            samplesConstant.push_back(proposedValue);
            densitiesConstant.push_back(proposedDensity);
            ++accRatio;
        } else {
            samplesConstant.push_back(samplesConstant.back());
            densitiesConstant.push_back(densitiesConstant.back());
        }
    }
}

template class RandomCalibration<RungeKuttaRandomH>;
template class RandomCalibration<RungeKuttaAddNoise>;