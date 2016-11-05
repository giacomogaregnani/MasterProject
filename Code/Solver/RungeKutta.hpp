#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

#include <Eigen/Dense>

using namespace Eigen;

// EULER FORWARD METHOD
class EulerForward {
private:
	int size;

	VectorXd (*f)(VectorXd, std::vector<double>&);

	std::vector<double> parameters;
public:
	EulerForward(int n, VectorXd (*func) (VectorXd, std::vector<double>&), std::vector<double> paramVec);

	VectorXd oneStep(VectorXd lastValue, double h);

	int getOrder(void);
};

// CLASSIC ORDER 4 RK METHOD
class RungeKutta {
private:
	int size;

	VectorXd (*f)(VectorXd, std::vector<double>&);

	std::vector<double> parameters;
public:	
	RungeKutta(int n, VectorXd (*func) (VectorXd, std::vector<double>&), std::vector<double> paramVec);

	VectorXd oneStep(VectorXd lastValue, double h);	

	int getOrder(void);
};


#endif
