#ifndef PARTICLEFILTER_HPP
#define PARTICLEFILTER_HPP

#include <EulerMaruyama.hpp>
#include <random>
#include <memory>
#include "ParFilTools.hpp"

class ParFil {
private:
    oneDimSde sde;
    double T;
    double IC;
    std::vector<double> obs;
    double noise;
    std::vector<std::vector<double>> X;
    std::vector<double> XOld;
    unsigned long nParticles;
    double eps;
    std::vector<double> W;
    double likelihood;
    std::default_random_engine seed;
    std::discrete_distribution<unsigned int> shuffler;
    std::vector<double> timeNoise;
    std::shared_ptr<ForwardPF> forwardSampler;

public:
    ParFil() = default;
    ~ParFil() = default;
    ParFil(std::vector<double>& y, double T, double IC,
           double noise, double eps, unsigned long nParticles,
           std::shared_ptr<ForwardPF>& forwardSampler);
    void compute(VectorXd& theta, std::vector<std::vector<double>>* mod = nullptr,
                 std::vector<double>* weights = nullptr);
    void computeDiffBridge(VectorXd& theta, std::vector<std::vector<double>>* mod = nullptr,
                           std::vector<double>* weights = nullptr, bool verbose = false);
    double getLikelihood() const;
    std::vector<double> sampleX();
    std::vector<std::vector<double>> getX() const;
};

#endif