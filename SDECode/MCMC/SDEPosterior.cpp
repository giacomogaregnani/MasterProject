#include <iostream>
#include <algorithm>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <chrono>
#include "SDEPosterior.hpp"
#include "../matplotlib-cpp-master/matplotlibcpp.h"

namespace plt = matplotlibcpp;

SDEPosterior::SDEPosterior(std::vector<double>& x, double T, double IC,
                           unsigned int sR, double noise, oneDimSde sde,
                           double eps, unsigned long M):
        T(T),
        IC(IC),
        samplingRatio(sR),
        obs(x),
        noise(noise),
        eps(eps),
        nMC(M)
{
    std::random_device dev;
    std::default_random_engine seed{dev()};
    Solver = EM1D(sde, seed);
    QOIobs = 0.0;
    for (auto it : x) {
        QOIobs += it;
    }
    QOIobs /= obs.size();
}

std::string GetLocalTime() {
    auto now(std::chrono::system_clock::now());
    auto seconds_since_epoch(
            std::chrono::duration_cast<std::chrono::seconds>(now.time_since_epoch()));

    // Construct time_t using 'seconds_since_epoch' rather than 'now' since it is
    // implementation-defined whether the value is rounded or truncated.
    std::time_t now_t(
            std::chrono::system_clock::to_time_t(
                    std::chrono::system_clock::time_point(seconds_since_epoch)));

    char temp[10];
    if (!std::strftime(temp, 10, "%H-%M-%S-", std::localtime(&now_t)))
        return "";

    return std::string(temp) +
           std::to_string((now.time_since_epoch() - seconds_since_epoch).count());
}

double SDEPosterior::computePosterior(VectorXd& theta)
{
    VectorXd thetaWithoutEps(theta.size()-1);
    thetaWithoutEps = theta.tail(theta.size()-1);
    double prior = -0.5 * thetaWithoutEps.dot(thetaWithoutEps);
    theta(0) = eps;

    double solution, likelihood;
    auto N = obs.size()-1;
    double h = T/N;
    auto nObs = static_cast<unsigned int>(std::round(N / samplingRatio));

    Solver.modifyParam(theta);

    /* std::vector<double> likelihoods(nMC);
    for (unsigned long i = 0; i < nMC; i++)
        likelihoods[i] = 0.0;

    std::vector<double> QOIVec(nMC);

    unsigned long k;
    #pragma omp parallel for num_threads(5) private(k)
    for (k = 0; k < nMC; k++) {
        solution = IC;
        for (unsigned long j = 0; j < nObs; j++) {
            for (unsigned long i = 0; i < samplingRatio; i++) {
                solution = Solver.oneStep(h, solution);
            }
            // dist = solution - obs[(j+1)*samplingRatio];
            // likelihoods[k] += -0.5 / (noise * noise) * dist * dist;
        }
        QOIVec[k] = solution;
    } */

    std::vector<double> solutions[N+1];
    std::vector<double> QOIVec(N+1);
    solution = IC;
    QOIVec[0] = IC;
    for (unsigned long j = 0; j < N; j++) {
            solution = Solver.oneStep(h, solution);
            QOIVec[j+1] = solution;
    }


    /* std::vector<double> timeVec(N+1);
    for (unsigned int i = 0; i < N+1; i++)
        timeVec[i] = h*i;
    plt::plot(timeVec, obs);
    plt::plot(timeVec, QOIVec);
    plt::show();

    plt::hist(obs, 20, "b", 0.3);
    plt::hist(QOIVec, 20, "r", 0.3);
    plt::show(); */

    // dist = solution - obs[(j+1)*samplingRatio];
    // likelihoods[k] += -0.5 / (noise * noise) * dist * dist;

    double minVal = std::min(*std::min_element(obs.begin(), obs.end()),
                             *std::min_element(QOIVec.begin(), QOIVec.end()));
    double maxVal = std::max(*std::max_element(obs.begin(), obs.end()),
                             *std::max_element(QOIVec.begin(), QOIVec.end()));

    unsigned long nF = 500;
    std::vector<double> FObs(nF), FX(nF);

    std::vector<double> sortObs = obs;
    std::sort(sortObs.begin(), sortObs.end());
    std::sort(QOIVec.begin(), QOIVec.end());

    std::vector<double> xValues(nF);
    double increment = (maxVal - minVal) / nF;
    for (unsigned int i = 0; i < nF; i++)
        xValues[i] = minVal + increment * i;

    auto itXOld = QOIVec.begin();
    auto itYOld = sortObs.begin();
    FX[0] = 0; FObs[0] = 0;
    std::vector<double>::iterator itXNew, itYNew;
    double distance = 0.0;
    for (unsigned int i = 1; i < nF; i++) {
        itXNew = std::upper_bound(itXOld, QOIVec.end(),  xValues[i]);
        itYNew = std::upper_bound(itYOld, sortObs.end(), xValues[i]);
        FX[i] = std::distance(itXOld, itXNew);
        FX[i] = FX[i-1] + FX[i] / QOIVec.size();
        FObs[i] = std::distance(itYOld, itYNew);
        FObs[i] = FObs[i-1] + FObs[i] / obs.size();
        itXOld = itXNew;
        itYOld = itYNew;
        distance = std::max(distance, std::abs(FX[i]-FObs[i]));
    }


    std::string str("/net/smana3/vol/vol2/anmc/garegnan/Desktop/Project/SDECode/plots/CDF");
    str = str + GetLocalTime();
    std::string plotTitle = std::to_string(distance) + " " + std::to_string(theta(1)) + " " + std::to_string(theta(2));
    plt::named_plot(plotTitle, xValues, FX);
    plt::plot(xValues, FObs);
    plt::legend();
    plt::save(str);
    plt::close();

    likelihood = -0.5 * distance * distance / (noise * noise);

    /* auto maxLikIt = std::max_element(likelihoods.begin(), likelihoods.end());
    double maxLik = *maxLikIt;
    likelihoods.erase(maxLikIt);

    double sum = 0;
    for (auto it : likelihoods) {
        sum += exp(it - maxLik);
    }
    likelihood = maxLik + std::log(1.0 + sum) - std::log(nMC); */

    return prior + likelihood;
}