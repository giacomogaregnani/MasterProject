#include "sRungeKutta.hpp"
#include <iostream>
#include <cmath>

detSROCK::detSROCK(odeDef ODE, std::vector<double> paramVec, int stages, double damping)
{
	f = ODE.odeFunc;
	size = ODE.size;
	parameters = paramVec;
	nStages = stages;
	omegaZero = 1.0 + damping / (stages * stages);
	chebCoeff.resize(nStages + 1);
	computeChebCoeff(omegaZero);
	stageCoeff.resize(nStages + 1);
	K.resize(nStages + 1);
	for (auto it : K) {
		it.resize(size);
	}
}

VectorXd detSROCK::oneStep(VectorXd solution, double h, int unUsed)
{
	// Hard-code the first two stages
	K[0] = stageCoeff[0][0] * solution;
	K[1] = stageCoeff[1][0] * solution + stageCoeff[1][1] * f(K[0], parameters);
	for (int i = 2; i < nStages + 1; i++) {
		K[i] = stageCoeff[i][0] * f(K[i - 1], parameters) + 
		       stageCoeff[i][1] * K[i - 1] + 
			   stageCoeff[i][2] * K[i - 2];
	}
	return K[nStages];
}

int detSROCK::getOrder(void)
{
	return 1;
}

void detSROCK::computeChebCoeff(double x)
{
	chebCoeff[0] = 1.0;
	chebCoeff[1] = x;
	for (int i = 2; i < nStages + 1; i++) {
		chebCoeff[i] = 2 * x * chebCoeff[i - 1] - chebCoeff[i - 2];
	}

	// Compute omegaOne (derivative of Cheb polynomials)
	std::vector<double> chebDerivatives(nStages + 1);
	chebDerivatives[0] = 0.0;
	chebDerivatives[1] = 1.0;
	chebDerivatives[2] = 4 * x;
	double doubI;
	for (int i = 3; i < nStages + 1; i++) {
		doubI = static_cast<double>(i);
		chebDerivatives[i] = doubI * (2 * chebCoeff[i - 1] + chebDerivatives[i - 2] / (doubI - 2));
	}
	omegaOne = chebCoeff.back() / chebDerivatives.back();
}

void detSROCK::computeStageCoeff(double h)
{
	// Hard-code the first two stages
	stageCoeff[0].resize(1);
	stageCoeff[0][0] = 1.0;
		
	stageCoeff[1].resize(2);
	stageCoeff[1][0] = 1.0;
	stageCoeff[1][1] = h * omegaOne / omegaZero;
	
	for (int i = 2; i < nStages + 1; i++) {
		stageCoeff[i].resize(3);
		stageCoeff[i][0] = 2 * h * omegaOne * chebCoeff[i - 1] / chebCoeff[i];
		stageCoeff[i][1] = 2 * omegaZero * chebCoeff[i - 1] / chebCoeff[i];
		stageCoeff[i][2] = -1.0 * chebCoeff[i - 2] / chebCoeff[i]; 
	}
}

// Classic first order RKC method
RKC::RKC(odeDef ODE, std::vector<double> paramVec,
		 int stages, double damping)
{
	f = ODE.odeFunc;
	size = ODE.size;
	parameters = paramVec;
	nStages = stages;
	kOld.resize(size);
	kNew.resize(size);
	kCurr.resize(size);
}

VectorXd RKC::oneStep(VectorXd lastValue, double h, int localStages)
{
	kOld = lastValue;
	double coeff = h / static_cast<double>(localStages * localStages);
	kNew = lastValue + f(lastValue, parameters) * coeff;
	coeff = 2.0 * coeff;

	for (int i = 2; i < localStages + 1; i++) {
		kCurr = f(kNew, parameters) * coeff +
			   kNew * 2.0 - kOld;
		kOld = kNew;
		kNew = kCurr;
	}

	return kCurr;
}

int RKC::getOrder(void)
{
	return 1;
}

void RKC::computeStageCoeff(double h) {}
