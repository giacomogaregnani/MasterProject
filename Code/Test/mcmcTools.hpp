#include <vector>
#include <Eigen/Dense>

using namespace Eigen;

// INSTRUMENTS (PRIOR AND LIKELIHOOD EVALUATIONS) (mcmcTools.cpp)
double evalPrior(VectorXd x, VectorXd mean, VectorXd variance, int size);

double evalLikelihood(std::vector<VectorXd> data,
					  std::vector<VectorXd> values,
					  int size, double variance);

double evalLogPrior(VectorXd x, VectorXd mean, VectorXd variance, int size);

double evalLogLikelihood(std::vector<VectorXd>& data,
						 std::vector<VectorXd>& values,
						 int size, double variance);

double evalSingleLikelihood(VectorXd data, VectorXd value, double variance);

// IMPLEMENTED IN THEIR CPP FILE
std::vector<VectorXd> metropolisHastings(VectorXd initialGuess, std::vector<double>& param, double sigma,
										 int size, double h, double finalTime, std::vector<VectorXd>& data,
										 std::vector<double>& dataTimes,
										 VectorXd (*odeFunc) (VectorXd, std::vector<double>&),
										 VectorXd priorMean, VectorXd priorVariance, int internalMC,
										 int nStepsMC, double varData, double* accRatio);

std::vector<VectorXd> sMetropolisHastings(VectorXd initialGuess, std::vector<double>& param, double sigma,
										  int size, double h, double finalTime,
										  std::vector<VectorXd>& data, std::vector<double>& dataTimes,
										  VectorXd (*odeFunc) (VectorXd, std::vector<double>&),
										  VectorXd priorMean, VectorXd priorVariance, int internalMC,
										  int nStepsMC, int nLev, double damping, double varData);

std::vector<VectorXd> MLmetropolisHastings(VectorXd initialGuess, std::vector<double>& param, double sigma,
										   int size, double h, double finalTime,
										   std::vector<VectorXd>& data,
										   std::vector<double>& dataTimes,
										   VectorXd (*odeFunc) (VectorXd, std::vector<double>&),
										   VectorXd priorMean, VectorXd priorVariance, int internalMC,
										   int nStepsMC, double varData, long int* costPerIteration,
										   double* accRatio);

std::vector<VectorXd> sMLmetropolisHastings(VectorXd initialCond, std::vector<double>& param, double sigma,
											int size, double h, double finalTime,
											std::vector<VectorXd>& data, std::vector<double>& dataTimes,
											VectorXd (*odeFunc) (VectorXd, std::vector<double>&),
											VectorXd priorMean, VectorXd priorVariance, int internalMC,
											int nStepsMC, double varData, double rho, double damping);

std::vector<VectorXd> gaussMetropolisHastings(VectorXd initialCond, std::vector<double>& param, double sigma,
											  int size, double h, double finalTime, std::vector<VectorXd>& data,
											  std::vector<double>& dataTimes,
											  VectorXd (*odeFunc) (VectorXd, std::vector<double>&),
											  VectorXd priorMean, VectorXd priorVariance,
											  int nStepsMC, double varData, double* accRatio);
