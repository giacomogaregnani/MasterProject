#include <cmath>
#include <vector>
#include <random>
#include <Solver.hpp>
#include <Eigen/Dense>
#include <ctime>
#include <iostream>
#include <fstream>
#include "problems.hpp"
#include "mcmcTools.hpp"
#include <Eigen/Eigenvalues>
#include <complex>
#include <chrono>

using namespace Eigen;

#define PI 3.1415926535897	

int main(int argc, char* argv[])
{
    // READ INPUT
	bool MC = false;
	bool stableMethod = false;
	bool multiLevel = false;

	if (argc > 1) {
		if (strcmp(argv[1], "-h") == 0 || argc < 5) {
			std::cout << "Usage: " << std::endl
					  << argv[0] << " namefile stable multi MC" << std::endl;
			return 0;
		}

		if (strcmp(argv[2], "stable") == 0) {
			stableMethod = true;
		}

		if (strcmp(argv[3], "multi") == 0) {
			multiLevel = true;
		}

		if (strcmp(argv[4], "mc") == 0) {
			MC = true;
		}
	}

    // TIME STEP AND NUMBER OF MC TRAJECTORIES
	std::vector<double> h = {};
	std::vector<int> nMCloop = {};

	if (MC && multiLevel) {
        int nExp = 4;
        double maxH = 0.125;
        for (int i = 0; i < nExp; i++) {
            h.push_back(maxH);
            maxH /= 2;
        }
        nMCloop.push_back(1);
	} else if (MC) {
		int nExp = 1;
		double maxH = 0.01;
        for (int i = 0; i < 100; i++) {
            nMCloop.push_back(650);
        }
        h.push_back(maxH);
        /* int nMC[3] = {300, 400, 1000};
        nMCloop.insert(nMCloop.end(), &nMC[0], &nMC[0] + 3);
		for (int i = 0; i < nExp; i++) {
			h.push_back(maxH);
			maxH /= 2;
        } */
	} else {
        int nExp = 4;
        double maxH = 0.001;
        for (int i = 0; i < nExp; i++) {
            h.push_back(maxH);
            maxH /= 2;
        }
		nMCloop.push_back(1);
	}

    // PROBLEM DATA
    odeDef odeModel;
    odeModel.ode = FITZNAG;
	setProblem(&odeModel);

	// Choose parameter list (on which inference will be done)
    std::vector<double> paramList = odeModel.refParam;

    // DATA ACQUISITION FROM REFSOL
    std::string refFile("refSolFitznag");
	std::fstream refSolution;
	std::string extension(".txt");
	std::string slash("/");
	std::string folder(DATA_PATH);
	refSolution.open(folder + slash + refFile + extension, std::ios_base::in);
	double finalTime;
	int nData;
	// read final time and number of data from reference solution
	refSolution >> finalTime;
	refSolution >> nData;
	std::vector<double> times(static_cast<unsigned int>(nData));
	std::vector<VectorXd> data(static_cast<unsigned int>(nData), VectorXd(odeModel.size));
	// read times where data was computed
	for (int i = 0; i < nData; i++) {
		refSolution >> times[i];
	}
	// read data
	for (int i = 0; i < nData; i++) {
		for (int j = 0; j < odeModel.size; j++) {
			refSolution >> data[i](j);
		}
	}

    // DATA UNCERTAINTY (COHERENT WITH REFSOL)
	double varData = 1e-2;

	time_t now;
	now = time(NULL);
	char charDate[18];
	if (now != -1) {
		strftime(charDate, 18, "%d_%m_%Y_%I_%M_", gmtime(&now));
	}
	std::string strDate(charDate);

    // a counter
    int count = 0;

    for (auto ith : h) {
		for (auto nInternalMC : nMCloop) {

			// PRIOR MEAN AND VARIANCE
			size_t nParam = paramList.size();	
			VectorXd priorMean(nParam), priorVariance(nParam);
			for (size_t i = 0; i < nParam; i++) {
				priorMean(i) = paramList[i];
				priorVariance(i) = 1.0;
			}
			
			// INITIAL MCMC GUESS
			std::vector<double> paramGuess(nParam);
			std::normal_distribution<double> disturbOnParam(0.0, 0.01);
			std::default_random_engine generator{(unsigned int) time(NULL)};
			for (size_t i = 0; i < nParam; i++) {
                paramGuess[i] = paramList[i];
				// paramGuess[i] = paramList[i] + disturbOnParam(generator);
			}

            // PARAMETERS OF THE CHAIN
			std::vector<VectorXd> mcmcPath;
			int nMCMC = 1000;
			double sigma = 0.5;

            // COST
            long int cost;
			double accRatio;

			StabValues stabParam;
			if (stableMethod) {
				stabParam.damping = 0.1;
                stabParam.method = stdRKC;
			}

            bool isPositive = false;
            if (odeModel.ode == BRUSS) {
                isPositive = true;
            }

            std::vector<double> likelihoods = {};

			if (MC) {
				if (!stableMethod) {
					if (!multiLevel) {
                        mcmcPath = MetropolisHastings(odeModel, paramGuess, sigma, ith, finalTime,
                                                      data, times, priorMean, priorVariance, nInternalMC,
                                                      nMCMC, varData, &cost, isPositive,
                                                      likelihoods, generator);
					} else {
						/* mcmcPath = MLmetropolisHastings(odeModel.initialCond, paramGuess, sigma, odeModel.size, ith, finalTime,
														data, times, odeModel.odeFunc, priorMean, priorVariance, nInternalMC,
														nMCMC, varData, &cost, &accRatio); */
					}
				} else {
					if (!multiLevel) {
						mcmcPath = sMetropolisHastings(odeModel, paramGuess, sigma, ith, finalTime,
													   data, times, priorMean, priorVariance, nInternalMC,
													   nMCMC, stabParam.damping, varData, &cost, isPositive,
                                                       likelihoods);
					} else {
						mcmcPath = sMLmetropolisHastings(odeModel, paramGuess, sigma, ith, finalTime,
														 data, times, priorMean, priorVariance,
														 nMCMC, varData, stabParam.damping, &cost);
					}
				}
			} else {
				mcmcPath = gaussMetropolisHastings(odeModel, paramGuess, sigma, ith, finalTime,
			                	                   data, times, priorMean, priorVariance,
							                       nMCMC, varData, &accRatio, stableMethod, stabParam);
			}

			// OUTPUT FILEs
            bool printResults = true;
            if (strcmp(argv[1], "debug") == 0) {
                printResults = false;
            }
            std::ofstream results;
			std::string hIter = std::to_string(static_cast<int>(ith * 1000));
            std::string nMC = std::to_string(nInternalMC);
            std::string nExperience = std::to_string(count++);
			std::string filepath;
			// default value to "test"
			if (argc > 1) {
				filepath = argv[1];
			} else {
				filepath = "test";
			}
			filepath = filepath + strDate + "_MC_" + nMC + "_h_" + hIter + "_Exp_" + nExperience + extension;
			std::string finalpath = folder + slash + filepath;
            if (printResults) results.open(finalpath , std::ofstream::out | std::ofstream::trunc);
			std::cout << finalpath << std::endl;

            // WRITE THE COST ON A FILE
            if (stableMethod && printResults && false) {
                std::ofstream costFile;
                filepath = argv[1] + strDate + "_MC_" + nMC + "_h_" + hIter + std::string("_cost") + extension;
                finalpath = folder + slash + filepath;
                costFile.open(finalpath, std::ofstream::out | std::ofstream::trunc);
                costFile << cost;
                costFile.close();
            }

			// WRITE RESULTS ON FILE
            if (printResults) {
                for (int i = 0; i < nMCMC; i++) {
                    results << (mcmcPath[i]).transpose() << " " << likelihoods[i] << "\n";
                }
            }

			std::cout << "acceptance Ratio = " << accRatio << std::endl;

            if (printResults) results.close();
		}
	}

	return 0;
}