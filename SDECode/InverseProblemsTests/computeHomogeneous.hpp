#ifndef COMPUTEHOMOGENEOUS_HPP
#define COMPUTEHOMOGENEOUS_HPP

#include <Eigen/Dense>
#include <vector>

using namespace Eigen;

std::vector<double> computeHomCoeffs(VectorXd&, double, double (*) (double));

std::vector<double> computeZs(double (*) (double), double, double);

#endif