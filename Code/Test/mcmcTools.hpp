#include <vector>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Solver.hpp>

using namespace Eigen;

// INSTRUMENTS (PRIOR AND LIKELIHOOD EVALUATIONS) (mcmcTools.cpp)
VectorXd computeMeans(std::vector<VectorXd>& input);

VectorXd computeVariances(std::vector<VectorXd>& input, VectorXd inputMean);

double phi(double x);

double evalPrior(VectorXd x, VectorXd mean, VectorXd variance, int size);

double evalLikelihood(std::vector<VectorXd> data,
					  std::vector<VectorXd> values,
					  int size, double variance);

double evalLogPrior(VectorXd x, VectorXd mean, VectorXd variance, int size);

double evalLogLikelihood(std::vector<VectorXd>& data,
						 std::vector<VectorXd>& values,
						 int size, double variance);

double evalSingleLikelihood(VectorXd data, VectorXd value, double variance);

double computeBatchMeans(std::vector<VectorXd>& data, double MCMean);

// IMPLEMENTED IN THEIR CPP FILE
std::vector<VectorXd> detMetropolisHastings(odeDef odeModel, std::vector<double>& param, double h,
											double finalTime, std::vector<VectorXd>& data, std::vector<double>& dataTimes,
											VectorXd priorMean, VectorXd priorVariance,
											int nStepsMC, double varData, long int* cost, bool posPar,
											std::vector<double>& likelihoods, std::default_random_engine& generator);

std::vector<VectorXd> MetropolisHastings(odeDef odeModel, std::vector<double>& param, double sigma,
										 double h, double finalTime, std::vector<VectorXd>& data, std::vector<double>& dataTimes,
										 VectorXd priorMean, VectorXd priorVariance, int internalMC,
										 int nStepsMC, double varData, long int* cost, bool posPar,
										 std::vector<double>& likelihoods, std::default_random_engine& generator);

void StoppedMetropolis(std::vector<VectorXd>& mcmcPath, odeDef odeModel, std::vector<double>& param, double sigma,
					   double h, double finalTime, std::vector<VectorXd>& data, std::vector<double>& dataTimes,
					   VectorXd priorMean, VectorXd priorVariance, int internalMC,
					   int nStepsMC, double varData, bool posPar,
					   std::vector<double>& likelihoods, std::default_random_engine& generator,
					   MatrixXd& S, int index);

void impStoppedMetropolis(std::vector<VectorXd>& mcmcPath,
                          odeDef odeModel, std::vector<double>& param, double sigma,
                          double h, double finalTime, std::vector<VectorXd>& data, std::vector<double>& dataTimes,
                          VectorXd priorMean, VectorXd priorVariance, int internalMC,
                          int nStepsMC, double varData, bool posPar,
                          std::vector<double>& likelihoods, std::default_random_engine& generator,
                          MatrixXd& S, int index, Butcher butcher);

std::vector<VectorXd> sMetropolisHastings(odeDef odeModel, std::vector<double>& param, double sigma,
										  double h, double finalTime, std::vector<VectorXd>& data, std::vector<double>& dataTimes,
										  VectorXd priorMean, VectorXd priorVariance, int internalMC,
										  int nStepsMC, double damping, double varData, long int* cost, bool posPar,
										  std::vector<double>& likelihoods);

/*std::vector<VectorXd> MLmetropolisHastings(VectorXd initialGuess, std::vector<double>& param, double sigma,
										   int size, double h, double finalTime,
										   std::vector<VectorXd>& data,
										   std::vector<double>& dataTimes,
										   VectorXd (*odeFunc) (VectorXd, std::vector<double>&),
										   VectorXd priorMean, VectorXd priorVariance, int internalMC,
										   int nStepsMC, double varData, long int* costPerIteration,
										   double* accRatio);*/

std::vector<VectorXd> sMLmetropolisHastings(odeDef odeModel, std::vector<double>& param, double sigma,
											double h, double finalTime,
											std::vector<VectorXd>& data,
											std::vector<double>& dataTimes,
											VectorXd priorMean, VectorXd priorVariance,
											int nStepsMC, double varData, double damping, long int* cost);

std::vector<VectorXd> gaussMetropolisHastings(odeDef odeModel, std::vector<double>& param, double sigma,
											  double h, double finalTime, std::vector<VectorXd>& data,
											  std::vector<double>& dataTimes,
											  VectorXd priorMean, VectorXd priorVariance, int nStepsMC,
											  double varData, double* accRatio, bool isStable, StabValues stab);

std::vector<VectorXd> testMetropolis(VectorXd oldGuess, int nIter, double* accRatio,
                                     double gamma, bool RAM, double desiredAlpha);

std::vector<VectorXd> testMetropolisConv(VectorXd oldGuess, int nIter, double* accRatio,
										 double gamma, bool RAM, double desiredAlpha,
										 MatrixXd& S, int* lastIndex, std::default_random_engine& generator);

MatrixXd RAMinit(double gamma, double desiredAlpha, int nParam);

MatrixXd RAMupdate(MatrixXd S, VectorXd w, double alpha, double desiredAlpha, int nParam, int nIt);

