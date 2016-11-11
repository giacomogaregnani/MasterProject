#include <cmath>
#include <iostream> 
#include "problems.hpp"
#include <Eigen/Core>

#ifndef PI
#define PI 3.1415926535897
#endif

// PROBLEMS
VectorXd lorenz(VectorXd argument, std::vector<double>& param) {
	double sigma = param[0], rho = param[1], beta = param[2];

	VectorXd result(3);
	result(0) = sigma * (argument(1) - argument(0));
	result(1) = argument(0) * (rho - argument(2)) - argument(1);
	result(2) = argument(0) * argument(1) - beta * argument(2);
	return result;
}

VectorXd test(VectorXd argument, std::vector<double>& param) {
	VectorXd result(3);
	for (int i = 0; i < 3; i++) {
		result(i) = param[i] * argument(i);
	}
	return result;
}

VectorXd fitznag(VectorXd argument, std::vector<double>& param)
{
	double a = param[0], b = param[1], c = param[2];
	VectorXd result(2);

	result(0) = c * (argument(0) - argument(0) * argument(0) * argument(0) / 3.0 + argument(1));
	result(1) = - 1.0 / c * (argument(0) - a + b * argument(1));

	return result;
}

MatrixXd fitznagJ(VectorXd argument, std::vector<double>& param)
{
    double b = param[1], c = param[2];
    MatrixXd result(2, 2);

    result(0, 0) = c * (1.0 - argument(0) * argument(0));
    result(0, 1) = c;
    result(1, 0) = -1.0 / c;
    result(1, 1) = -1.0 * b / c;

    return result;
}

VectorXd testOneD(VectorXd argument, std::vector<double>& param) {
	return argument * param[0];
}

MatrixXd testOneDJ(VectorXd argument, std::vector<double>& param) {
	return MatrixXd::Identity(1, 1) * param[0];
}

VectorXd vdPol(VectorXd argument, std::vector<double>& param) {
	VectorXd result(2);
	result(0) = argument(1);
	result(1) = param[0] * (1 - argument(0) * argument(0)) * argument(1) - argument(0);
	return result;
}


MatrixXd vdPolJ(VectorXd argument, std::vector<double>& param) {
	MatrixXd result(2, 2);
	result(0, 0) = 0.0;
	result(0, 1) = 1.0;
	result(1, 0) = -2.0 * param[0] * argument(1) * argument(0) - 1.0;
	result(0, 1) = param[0] * (1 - argument(0) * argument(0));
	return result;
}

VectorXd robertson(VectorXd argument, std::vector<double>& param) {
	VectorXd result(3);
	result(0) = -0.04 * argument(0) + 10000 * argument(1) * argument(2);
	result(2) = 0.04 * argument(0) - 10000 * argument(1) * argument(3) - 3 * 10e7 * argument(3) * argument(3);
	result(3) = 3 * 10e7 * argument(3) * argument(3);

	return result;
}

VectorXd bruss(VectorXd argument, std::vector<double>& param)
{
	double alpha = param[0] / 50.0;
	int N = static_cast<int>(argument.size());
	int Nh = N / 2;
	double dX = 1.0 / (Nh + 1);
	double coeff = alpha / (dX * dX);
	VectorXd result(N);
	
	// u species
	result(0) = 1.0 + argument(0) * argument(0) * argument(Nh) - 4.0 * argument(0)
	          + coeff * (1.0 - 2 * argument(0) + argument(1));
	for (int i = 1; i < Nh - 1; i++) {
		result(i) = 1.0 + argument(i) * argument(i) * argument(Nh + i) 
			  - 4.0 * argument(i) + coeff * (argument(i - 1) - 2.0 * argument(i) + argument(i + 1));			
	}
	result(Nh - 1) = 1.0 + argument(Nh - 1) * argument(Nh - 1) * argument(N - 1) - 4.0 * argument(Nh - 1)
	          + coeff * (argument(Nh - 2) - 2.0 * argument(Nh - 1) + 1.0); 

	// v species
	result(Nh) = 3.0 * argument(0) - argument(0) * argument(0) * argument(Nh)
	          + coeff * (3.0 - 2.0 * argument(Nh) + argument(Nh + 1)); 
	for (int i = Nh + 1; i < N - 1; i++) {
		result(i) = 3.0 * argument(i - Nh) - argument(i - Nh) * argument(i - Nh) * argument(i)
			  + coeff * (argument(i - 1) - 2 * argument(i) + argument(i + 1)); 
	}
	result(N - 1) = 3.0 * argument(Nh - 1) - argument(Nh - 1) * argument(Nh - 1) * argument(N - 1)
	          + coeff * (argument(N - 2) - 2 * argument(N - 1) + 3.0); 

	return result;	
}

MatrixXd brussJ(VectorXd argument, std::vector<double>& param)
{
	double alpha = param[0] / 50.0;
	int N = static_cast<int>(argument.size());
	int Nh = N / 2;
	double dX = 1.0 / (Nh + 1);
	double coeff = alpha / (dX * dX);

	// Initialize the matrix
	MatrixXd result = MatrixXd::Zero(N, N);

	// du'_i / du_j and du'_i / dv_j
	result(0, 0) = 2.0 * argument(0) * argument(Nh) - 4.0 - 2.0 * coeff;
	result(0, 1) = coeff;
	result(0, Nh) = argument(0) * argument(0);
	for (int i = 1; i < Nh - 1; i++) {
		result(i, i) = 2.0 * argument(i) * argument(Nh + i) - 4.0 - 2.0 * coeff;
		result(i, i - 1) = coeff;
		result(i, i + 1) = coeff;
		result(i, i + Nh) = argument(i) * argument(i);
	}
	result(Nh - 1, Nh - 1) = 2.0 * argument(Nh - 1) * argument(N - 1) - 4.0 - 2.0 * coeff;
	result(Nh - 1, Nh - 2) = coeff;
	result(Nh - 1, N - 1) = argument(Nh - 1) * argument(Nh - 1);

	// dv'_i / du_j and dv'_i / dv_j
	result(Nh, 0) = 3.0 - 2.0 * argument(0) * argument(Nh);
	result(Nh, Nh) = -1.0 * argument(0) * argument(0) - 2.0 * coeff;
	result(Nh, Nh + 1) = coeff;
	for (int i = Nh + 1; i < N - 1; i++) {
		result(i, i) = -1.0 * argument(i - Nh) * argument(i - Nh) - 2.0 * coeff;
		result(i, i - 1) = coeff;
		result(i, i + 1) = coeff;
		result(i, i - Nh) = 3.0 - 2 * argument(i - Nh) * argument(i);
	}
	result(N - 1, N - 1) = -1.0 * argument(Nh - 1) * argument(Nh - 1) - 2 * coeff;
	result(N - 1, N - 2) = coeff;
	result(N - 1, Nh - 1) = 3.0 - 2 * argument(Nh - 1) * argument(N - 1);

	return result;
}

// For Poisson test, a class avoids to form the matrix only once
MatrixXd Poisson::A(10, 10);

VectorXd Poisson::poissonTest(VectorXd argument, std::vector<double>& param) {
	long int N = argument.size();
	VectorXd result(N);
    return A * argument * param[0];
}

Poisson::Poisson(int N) {
    A = MatrixXd::Zero(N, N);
    for (int i = 0; i < N - 1; i++) {
        A(i, i) = -2.0;
        A(i, i + 1) = 1.0;
        A(i + 1, i) = 1.0;
    }
    A(N - 1, N - 1) = -2.0;
}

void setProblem(odeDef* odeModel)
{
	std::vector<double> a = {1.0, 2.0, 3.0};
	switch ((*odeModel).ode) {
		case FITZNAG:
			(*odeModel).size = 2;
			((*odeModel).initialCond).resize((*odeModel).size);
			((*odeModel).initialCond)(0) = -1.0;
			((*odeModel).initialCond)(1) = 1.0;
			(*odeModel).odeFunc = &fitznag;
            (*odeModel).odeJac = &fitznagJ;
			break;
		case LORENZ:
			((*odeModel).size) = 3;
			((*odeModel).initialCond).resize((*odeModel).size);
			((*odeModel).initialCond)(0) = -10.0;
			((*odeModel).initialCond)(1) = -1.0;
			((*odeModel).initialCond)(2) = 40.0;
			(*odeModel).odeFunc = &lorenz;
			break;	
		case TEST:
			(*odeModel).size = 3;
			((*odeModel).initialCond).resize((*odeModel).size);
			((*odeModel).initialCond)(0) = 1.0;
			((*odeModel).initialCond)(1) = 2.0;
			((*odeModel).initialCond)(2) = 3.0;
			(*odeModel).odeFunc = &test;
			break;
		case TEST1D:
			(*odeModel).size = 1; 
			((*odeModel).initialCond).resize((*odeModel).size);
			((*odeModel).initialCond)(0) = 1.0;
			(*odeModel).odeFunc = &testOneD;
			(*odeModel).odeJac = &testOneDJ;
			break;
		case VDPOL:
			(*odeModel).size = 2;
			((*odeModel).initialCond).resize((*odeModel).size);
			((*odeModel).initialCond)(0) = 2.0;
			((*odeModel).initialCond)(1) = 0.0; 
			(*odeModel).odeFunc = &vdPol;
			(*odeModel).odeJac = &vdPolJ;
			break;
		case ROBERTSON:
			(*odeModel).size = 3;
			((*odeModel).initialCond).resize((*odeModel).size);
			((*odeModel).initialCond)(0) = 1.0;
			((*odeModel).initialCond)(1) = 0.0;
			((*odeModel).initialCond)(2) = 0.0;
			(*odeModel).odeFunc = &robertson;
			break;
		case BRUSS:
			(*odeModel).size = 10;
			((*odeModel).initialCond).resize((*odeModel).size);
			double x;
			for (int i = 0; i < (*odeModel).size / 2; i++) {
				x = static_cast<double>(i + 1) / ((*odeModel).size / 2.0 + 1);
				((*odeModel).initialCond)(i) = 1.0 + 0.5 * sin(2 * PI * x);	
			}
			for (int i = (*odeModel).size / 2; i < (*odeModel).size; i++) {
				((*odeModel).initialCond)(i) = 3.0;
			}
			(*odeModel).odeFunc = &bruss;
			(*odeModel).odeJac = &brussJ;
			break;
		case POISSON:
            Poisson P(10);
			(*odeModel).size = 10;
			((*odeModel).initialCond).resize((*odeModel).size);
			for (int i = 0; i < (*odeModel).size; i++) {
				((*odeModel).initialCond)(i) = 1.0;
			}
			(*odeModel).odeFunc = &P.poissonTest;
			break;
	}
}


