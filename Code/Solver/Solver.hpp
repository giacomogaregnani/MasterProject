#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <cmath>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <random>
#include <ctime>
#include "RungeKutta.hpp"
#include "sRungeKutta.hpp"
#include "Structures.hpp"
#include <unsupported/Eigen/KroneckerProduct>

// SUPPORT STRUCTURES AND TOOLS (SEE PROBLEMS.HPP)

double powerMethod(MatrixXd A, double tol, int nMax);

MatrixXd triInv(MatrixXd A, int size);

// PROBABILISTIC METHOD 
template <class T> 
class ProbMethod { 
private:
	std::shared_ptr<T> detSolver;

	int size;

	VectorXd solution;

	double sigma;

	double rootsigma;	

	std::normal_distribution<double> normalDist;	

	double h;

	double hfunc;

	VectorXd IC;

	std::vector<double> parameters;
public:
	ProbMethod(odeDef ODE, double step, std::vector<double> paramVec,
			   double stoch);

	VectorXd& getSolution(void);

	void oneStep(std::default_random_engine& generator, double step);

	VectorXd oneStepGiven(VectorXd& xi, double step);

	void resetIC(void);
};

class impProbMethod {
private:
	std::shared_ptr<ImplicitRK> detSolver;

	int size;

	VectorXd solution;

	double sigma;

	double rootsigma;

	std::normal_distribution<double> normalDist;

	double h;

	double hfunc;

	VectorXd IC;

	std::vector<double> parameters;
public:
	impProbMethod(odeDef ODE, double step, std::vector<double> paramVec,
			   double stoch, MatrixXd A, VectorXd b, int order);

	VectorXd& getSolution(void);

	void oneStep(std::default_random_engine& generator, double step);

	VectorXd oneStepGiven(VectorXd& xi, double step);

	void resetIC(void);
};


// STABILIZED PROBABILISTIC METHOD 
template <class T> 
class sProbMethod {
private:
	std::shared_ptr<T> detSolver;

	int size;

	VectorXd solution;

	double sigma;

	double rootsigma;	

	std::normal_distribution<double> normalDist;	

	double h;

	double hfunc;

	VectorXd IC;

	std::vector<double> parameters;

	int nStages;
public:
	sProbMethod(int n, double step,  
		   VectorXd initialCond, std::vector<double> paramVec,
		   VectorXd (*func) (VectorXd, std::vector<double>&),
		   double stoch, int nRCKStages, double damping);

	VectorXd& getSolution(void);

	void oneStep(std::default_random_engine& generator, double step, int localStages);

	VectorXd oneStepGiven(VectorXd& xi, double step, int localStages);

	void resetIC(void);
};

// MULTI-LEVEL MONTE CARLO
/*template <class T>
 class MLMC {
private:
	// Solver
	std::vector<std::shared_ptr<ProbMethod<T>>> solver;

	int detOrder;

	// Problem size
	int size;

	// Final time solution
	std::vector<VectorXd> solution;
	
	// Initial condition
	VectorXd initialCond;

	// Wanted accuracy
	double epsilon;	

	// Number of levels
	int nLevels;

	// Number of trajectories
	std::vector<long int> nTrajectories;

	// Timesteps 
	std::vector<double> timesteps;	
	std::vector<int> nSteps;
			
	// Final Time
	double finalTime;

	// generators
	std::normal_distribution<double> normalDist;
	std::default_random_engine generator{(unsigned int) time(NULL)};

	// function Phi
	double (*phiMLMC) (VectorXd, VectorXd, double);

	// only for MCMC application
	VectorXd dataMCMC;
	double varDataMCMC;

public:
	MLMC(int n, double accuracy, VectorXd initialCondition, std::vector<double> parameters,
		 VectorXd (*func) (VectorXd, std::vector<double>&), double stoch,
	     int order, double finalTime, bool fromTimestep, double hL,
	     double (*phi) (VectorXd, VectorXd, double),
         VectorXd data, double varianceData);

	long int cost(void);	

	double compute(void);

	double gethL(void);
};*/

// STABILIZED MULTI-LEVEL MONTE CARLO
template <class T>
class sMLMC {
private:
	// Solver
	std::vector<std::shared_ptr<sProbMethod<T>>> solver;
	int detOrder;

	// number of stages
	std::vector<int> nStages;

	// Problem size
	int size;

	// Final time solution
	std::vector<VectorXd> solution;
	
	// Initial condition
	VectorXd initialCond;

	// Wanted accuracy
	double epsilon;	

	// Number of levels
	int nLevels;

	// Number of trajectories
	std::vector<long int> nTrajectories;

	// Timesteps 
	std::vector<double> timesteps;	
	std::vector<int> nSteps;
			
	// Final Time
	double finalTime;

	// generators
	std::normal_distribution<double> normalDist;	
	std::default_random_engine generator{(unsigned int) time(NULL)};
	
	// function Phi
	double (*phiMLMC) (VectorXd, VectorXd, double);

	// only for MCMC application
	VectorXd dataMC;
	double varData;

public:
	sMLMC(int n, double accuracy, VectorXd initialCondition, std::vector<double> parameters,
          VectorXd (*func) (VectorXd, std::vector<double>&), double stoch,
          int order, double finalTime, double damping, double probStiffness,
          bool fromTimestep, double hL,
          double (*phi) (VectorXd, VectorXd, double),
          VectorXd data, double varianceData);

	long int cost(void);	

	double compute(void);

	double gethL(void);
};

// Third order Gauss cubature method
class ThirdOrderGauss {
	// system size
	int size;

	// parameters defining the system
	std::vector<double> theta;

	// The ode
	odeDef ODE;

	// Variance and square root of the variance
    MatrixXd P;
	MatrixXd lChol;
	MatrixXd sqrtP;
	void updateSqrtP(void);
	MatrixXd xAntiSym;

    // Correct the fact that sqrt(zeros) does not work in Eigen (-nan -nan ...)
    MatrixXd matrixSqrt(MatrixXd&);

	// Mean and old value of the mean
    VectorXd m;
	VectorXd oldM;

	// drift function
	VectorXd (*f) (VectorXd, std::vector<double>&);
    MatrixXd (*Fx)(VectorXd, std::vector<double>&);

	// diffusion function
	MatrixXd (*L) (VectorXd, std::vector<double>&, double, double);
	double h;
	double sigma;

	// sigma points
	std::vector<VectorXd> xi;
	std::vector<VectorXd> computeSigmaPoints(void);

	// data variance
	MatrixXd varData;

	// private EFupdates
	VectorXd meanUpdateFct(VectorXd, MatrixXd);
	MatrixXd varianceUpdateFct(VectorXd, MatrixXd);
    VectorXd TmeanUpdateFct(VectorXd);
    MatrixXd TvarianceUpdateFct(VectorXd, MatrixXd);
	MatrixXd TsqrtVarianceUpdateFct(VectorXd mIntern, MatrixXd LIntern);
	void EFupdates(int nSteps);
    void stabUpdates(int nSteps);
	void KalmanUpdate(VectorXd data);

	// evaluate Gaussian
	double evaluateGaussian(VectorXd data);

	// stability parameters
    StabValues stabParam;
	bool stableMethod;

    // Stability tools
    void computeStageCoeff(void);
    void computeChebCoeff(double x);
    double omegaZero;
    double omegaOne;
    std::vector<double> chebCoeff;
    std::vector<std::vector<double>> stageCoeff;
    std::vector<VectorXd> kMean;
    std::vector<MatrixXd> kVar;

public:
	ThirdOrderGauss(odeDef odeModel, MatrixXd initialVariance,
					std::vector<double> param,
					MatrixXd (*diffusion) (VectorXd, std::vector<double>&, double, double),
					double dataVariance, double step, double sigmadiff, bool isStiff,
					StabValues* stability);

	VectorXd getMean(void);
	MatrixXd getVariance(void);

	double oneStep(VectorXd data, int nSteps);
};

#endif