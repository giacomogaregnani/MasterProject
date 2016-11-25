#ifndef SRUNGEKUTTA_H
#define SRUNGEKUTTA_H

#include <Eigen/Dense>

using namespace Eigen;

// DETERMINISTIC PART OF THE S-ROCK METHOD
class detSROCK {
private:
	int size;

	VectorXd (*f)(VectorXd, std::vector<double>&);

	std::vector<double> parameters;

	int nStages;

	std::vector<double> chebCoeff;

	std::vector<std::vector<double>> stageCoeff;

	double omegaZero;

	double omegaOne;

	std::vector<VectorXd> K;
public:	
	detSROCK(int n, VectorXd (*func) (VectorXd, std::vector<double>&), 
                 std::vector<double> paramVec, int nRKCStages, double damping);

	VectorXd oneStep(VectorXd lastValue, double h, int localStages);

	int getOrder(void);

	void computeChebCoeff(double x);

	void computeStageCoeff(double h);
};

class RKC {
private:
	int size;

	VectorXd (*f)(VectorXd, std::vector<double>&);

	std::vector<double> parameters;

	int nStages;

	VectorXd kOld;
	VectorXd kNew;
	VectorXd kCurr;

public:
	RKC(int n, VectorXd (*func) (VectorXd, std::vector<double>&),
			 std::vector<double> paramVec, int nRKCStages, double damping);

	VectorXd oneStep(VectorXd lastValue, double h, int localStages);

	int getOrder(void);

	void computeStageCoeff(double h);
};


#endif
