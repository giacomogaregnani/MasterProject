#ifndef READOBSERVATIONS_HPP
#define READOBSERVATIONS_HPP

#include <Eigen/Dense>
#include <fstream>
#include <vector>

using namespace Eigen;

void ReadObservations(std::vector<double>& tObs,
                      std::vector<VectorXd>& observations,
                      unsigned int nObs, unsigned int size,
                      std::string& filename);

#endif