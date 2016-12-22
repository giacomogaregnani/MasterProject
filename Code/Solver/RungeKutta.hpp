#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

#include <Eigen/Dense>
#include "Structures.hpp"

using namespace Eigen;

VectorXd newton(odeDef ode, VectorXd uOld, std::vector<double>& param,
                double h, int maxIter, double tol);

VectorXd newtonGeneral(odeDef ode, VectorXd uOld, std::vector<double>& param,
					   double h, int maxIter, double tol, MatrixXd A, VectorXd b);

// EULER FORWARD METHOD
class EulerForward {
private:
	int size;

	VectorXd (*f)(VectorXd, std::vector<double>&);

	std::vector<double> parameters;

    odeDef theODE;
public:
	EulerForward(odeDef ODE, std::vector<double> paramVec);

	VectorXd oneStep(VectorXd lastValue, double h);

	int getOrder(void);
};

// EULER BACKWARDS METHOD
class EulerBackwards {
private:
	int size;

	VectorXd (*f)(VectorXd, std::vector<double>&);

	std::vector<double> parameters;

    odeDef theODE;
public:
	EulerBackwards(odeDef ODE, std::vector<double> paramVec);

	VectorXd oneStep(VectorXd lastValue, double h);

	int getOrder(void);
};

// CLASSIC ORDER 4 RK METHOD
class RungeKutta {
private:
	int size;

	VectorXd (*f)(VectorXd, std::vector<double>&);

	std::vector<double> parameters;

    odeDef theODE;
public:	
	RungeKutta(odeDef ODE, std::vector<double> paramVec);

	VectorXd oneStep(VectorXd lastValue, double h);	

	int getOrder(void);
};

// MIDPOINT METHOD
class MidPoint {
private:
	int size;

	VectorXd (*f)(VectorXd, std::vector<double>&);

	std::vector<double> parameters;

	odeDef theODE;
public:
	MidPoint(odeDef ODE, std::vector<double> paramVec);

	VectorXd oneStep(VectorXd lastValue, double h);

	int getOrder(void);
};

// GENERAL IMPLICIT RUNGE KUTTA METHOD
class ImplicitRK {
private:
	odeDef ODE;

	std::vector<double> parameters;

	MatrixXd A;

	VectorXd b;

	int nStages;

	int order;
public:
	ImplicitRK(odeDef inODE, std::vector<double> paramVec,
			   MatrixXd inA, VectorXd inB, int inOrder);

	VectorXd oneStep(VectorXd lastValue, double h);

	int getOrder(void);
};


#endif
