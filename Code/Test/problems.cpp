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

VectorXd fitznag(VectorXd argument, std::vector<double>& param) {
	double a = param[0], b = param[1], c = param[2];
	
	VectorXd result(2);
	result(0) = c * (argument(0) - argument(0) * argument(0) * argument(0) / 3.0 + argument(1));
	result(1) = - 1.0 / c * (argument(0) - a + b * argument(1));

	return result;
}

VectorXd testOneD(VectorXd argument, std::vector<double>& param) {
	return param[0] * argument;
}

VectorXd vdPol(VectorXd argument, std::vector<double>& param) {
	VectorXd result(2);
	result(0) = argument(1);
	result(1) = param[0] * (1 - argument(0) * argument(0)) * argument(1) - argument(0);

	return result;	
}

VectorXd robertson(VectorXd argument, std::vector<double>& param) {
	VectorXd result(3);
	result(0) = -0.04 * argument(0) + 10000 * argument(1) * argument(2);
	result(2) = 0.04 * argument(0) - 10000 * argument(1) * argument(3) - 3 * 10e7 * argument(3) * argument(3);
	result(3) = 3 * 10e7 * argument(3) * argument(3);

	return result;
}

VectorXd bruss(VectorXd argument, std::vector<double>& param) {
	double alpha = param[0];
	double dX = param[1];
	double coeff = alpha / (dX * dX);
	int N = static_cast<int>(argument.size());
	int Nh = N / 2;
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

// For Poisson test, a class avoids to form the matrix only once
MatrixXd Poisson::A(10, 10);

VectorXd Poisson::poissonTest(VectorXd argument, std::vector<double>& param) {
	long int N = argument.size();
	VectorXd result(N);
    // std::cout << A << std::endl;
    return param[0] * A * argument;
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


void setProblem(VectorXd* initialCond, VectorXd (**odeFunc) (VectorXd, std::vector<double>&), 
	        problems choice, int* size)
{
	std::vector<double> a = {1.0, 2.0, 3.0};
	switch (choice) {
		case FITZNAG:
			*size = 2;
			(*initialCond).resize(*size);
			(*initialCond)(0) = -1.0;
			(*initialCond)(1) = 1.0;
			*odeFunc = &fitznag;
			break;
		case LORENZ:
			(*size) = 3;
			(*initialCond).resize(*size);
			(*initialCond)(0) = -10.0;
			(*initialCond)(1) = -1.0;
			(*initialCond)(2) = 40.0;
			*odeFunc = &lorenz;
			break;	
		case TEST:
			*size = 3;
			(*initialCond).resize(*size);
			(*initialCond)(0) = 1.0;
			(*initialCond)(1) = 2.0;
			(*initialCond)(2) = 3.0;
			*odeFunc = &test;
			//std::cout << odeFunc(*initialCond, a) << std::endl;
			std::cout << odeFunc << std::endl;
			break;
		case TEST1D:
			*size = 1; 
			(*initialCond).resize(*size);
			(*initialCond)(0) = 1.0;
			*odeFunc = &testOneD;
			break;
		case VDPOL:
			*size = 2;
			(*initialCond).resize(*size);
			(*initialCond)(0) = 2.0;
			(*initialCond)(1) = 0.0; 
			*odeFunc = &vdPol;
			break;
		case ROBERTSON:
			*size = 3;
			(*initialCond).resize(*size);
			(*initialCond)(0) = 1.0;
			(*initialCond)(1) = 0.0;
			(*initialCond)(2) = 0.0;
			*odeFunc = &robertson;
			break;
		case BRUSS:
			*size = 100;
			(*initialCond).resize(*size);
			double x;
			for (int i = 0; i < *size / 2; i++) {
				x = static_cast<double>(i + 1) / static_cast<double>(*size / 2.0 + 1);
				(*initialCond)(i) = 1.0 + 0.5 * sin(2 * PI * x);	
			}
			for (int i = *size / 2; i < *size; i++) {
				(*initialCond)(i) = 3.0;
			}
			*odeFunc = &bruss;			
			break;
		case POISSON:
            Poisson P(10);
			*size = 10;
			(*initialCond).resize(*size);
			for (int i = 0; i < *size; i++) {
				(*initialCond)(i) = 1.0;
			}
			*odeFunc = &P.poissonTest;
			break;
	}
}


