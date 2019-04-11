#ifndef GENERATEOBSERVATIONS_HPP
#define GENERATEOBSERVATIONS_HPP

#include <EulerMaruyama.hpp>

std::vector<double> generateObservations1D(oneDimSde sde, double IC, VectorXd& param,
                                           double T, unsigned int N);

#endif