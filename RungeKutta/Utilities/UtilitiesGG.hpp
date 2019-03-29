#ifndef UTILITIESGG_H
#define UTILITIESGG_H

#include <iostream>
#include <fstream>
#include <string>
#include <memory>

#include <Eigen/Dense>

using namespace Eigen;

VectorXd loadRefSolution(std::string& fileName, int solSize);

void loadObservations(std::vector<double>& obsTimes, std::vector<VectorXd>& obs,
                      std::string& fileName, int solSize, double* noise);

void loadObservationsHeat(VectorXd& obs, std::string& fileName, int* solSize, double* noise);

void computeOrderOfConvergence(std::vector<double>& errors,
                               std::vector<double>& orders,
                               double ratio);

void printConvergenceInfo(std::vector<double>& errors,
                          std::vector<double>& orders);

VectorXd normalZeroMeanRandVec(int size,
                               std::default_random_engine* generator,
                               double stdDev);

std::vector<double> EigVecToStdVec(VectorXd& vec);

VectorXd StdVecToEig(std::vector<double>& vec);

#endif //UTILITIESGG_H