#ifndef PARTICLEFILTERFULL_HPP
#define PARTICLEFILTERFULL_HPP

#include <EulerMaruyama.hpp>
#include <random>
#include <memory>
#include "ParFilTools.hpp"

class ParFilMod {
private:
    double T;
    double IC;
    std::vector<double> obs;
    double noise;
    std::vector<std::vector<Vector2d>> X;
    std::vector<Vector2d> XOld;
    unsigned long nParticles;
    double eps;
    std::vector<double> W;
    double likelihood;
    std::default_random_engine seed;
    std::discrete_distribution<unsigned int> shuffler;
    std::shared_ptr<ForwardPFModErr> forwardSampler;

public:
    ParFilMod() = default;
    ~ParFilMod() = default;
    ParFilMod(std::vector<double> &y, double T, double IC,
               double noise, double eps, unsigned long nParticles,
               std::shared_ptr<ForwardPFModErr>& forwardSampler);
    void compute(VectorXd& theta, bool verbose = false);
    void computeIS(VectorXd& theta, bool verbose = false);
    double getLikelihood() const;
    std::vector<Vector2d> sampleX();
    std::vector<std::vector<Vector2d>> getX() const;
    std::vector<double> getW() const;
};

#endif