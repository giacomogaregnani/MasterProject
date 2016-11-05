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
        int nExp = 6;
        double maxH = 0.1;
        for (int i = 0; i < nExp; i++) {
            h.push_back(maxH);
            maxH /= 2;
        }
        nMCloop.push_back(1);
	} else if (MC) {
        int nExp = 6;
        double maxH = 0.05;
        for (int i = 0; i < nExp; i++) {
            h.push_back(maxH);
            maxH /= 2;
        }
        nMCloop.push_back(100);
        // nMCloop.push_back(1000);
		//int nDouble[] = {1, 10, 100, 1000, 10000};
		//nMCloop.insert(nMCloop.end(), nDouble, nDouble + 5);
	} else {
        int nExp = 6;
        double maxH = 0.1;
        for (int i = 0; i < nExp; i++) {
            h.push_back(maxH);
            maxH /= 2;
        }
		nMCloop.push_back(1);
	}

    // PROBLEM DATA
	problems problem = LORENZ;
	VectorXd (*odeFunc) (VectorXd, std::vector<double>&);
	int size;
	VectorXd initialCond;
	setProblem(&initialCond, &odeFunc, problem, &size);

	// Choose parameter list (on which inference will be done)
    std::vector<double> paramList = {10.0, 28.0, 8.0 / 3.0};
	// std::vector<double> paramList = {1.0};

    // DATA ACQUISITION FROM REFSOL
    std::string refFile("refSolLorenz.txt");
	std::fstream refSolution;
	std::string extension(".txt");
	std::string slash("/");
	std::string folder(DATA_PATH);
	refSolution.open(folder + slash + refFile, std::ios_base::in);
	double finalTime;
	int nData;
	// read final time and number of data from reference solution
	refSolution >> finalTime;
	refSolution >> nData;
	std::vector<double> times(static_cast<unsigned int>(nData));
	std::vector<VectorXd> data(static_cast<unsigned int>(nData), VectorXd(size));
	// read times where data was computed
	for (int i = 0; i < nData; i++) {
		refSolution >> times[i];
	}
	// read data
	for (int i = 0; i < nData; i++) {
		for (int j = 0; j < size; j++) {
			refSolution >> data[i](j);
		}
	}

    // DATA UNCERTAINTY (COHERENT WITH REFSOL)
	double varData = 1e-3;

    // Counter and time for output files
    int count = 0;

	time_t now;
	now = time(NULL);
	char charDate[18];
	if (now != -1) {
		strftime(charDate, 18, "%d_%m_%Y_%I_%M", gmtime(&now));
	}
	std::string strDate(charDate);

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
				paramGuess[i] = paramList[i] + 0.1 * disturbOnParam(generator);
			}

            // PARAMETERS OF THE CHAIN
			std::vector<VectorXd> mcmcPath;
			int nMCMC = 100000;
			double sigma = 0.5;

            // ONLY FOR STABLE METHODS
			double damping = 0.0;
			double rho = 25;
			int nStages = static_cast<int>(pow(3.0 * rho * ith, 0.5));

            // COST
            long int cost;
			double accRatio;

			if (MC) {
				if (!stableMethod) {
					if (!multiLevel) {
						mcmcPath = metropolisHastings(initialCond, paramGuess, sigma, size, ith, finalTime, 
													  data, times, odeFunc, priorMean, priorVariance, nInternalMC,
													  nMCMC, varData, &accRatio);
					} else {
						mcmcPath = MLmetropolisHastings(initialCond, paramGuess, sigma, size, ith, finalTime, 
														data, times, odeFunc, priorMean, priorVariance, nInternalMC,
														nMCMC, varData, &cost, &accRatio);
					}
				} else {
					if (!multiLevel) {
						mcmcPath = sMetropolisHastings(initialCond, paramGuess, sigma, size, ith, finalTime, 
													   data, times, odeFunc, priorMean, priorVariance, nInternalMC,
													   nMCMC, nStages, damping, varData);
					} else {
						mcmcPath = sMLmetropolisHastings(initialCond, paramGuess, sigma, size, ith, finalTime, 
														 data, times, odeFunc, priorMean, priorVariance, nInternalMC,
														 nMCMC, varData, rho, damping);
					}
				}
			} else {
				mcmcPath = gaussMetropolisHastings(initialCond, paramGuess, sigma, size, ith, finalTime,
			                	                   data, times, odeFunc, priorMean, priorVariance, 
							                       nMCMC, varData, &accRatio);
			}

			// OUTPUT FILEs
            bool printResults = true;
            if (strcmp(argv[1], "debug") == 0) {
                printResults = false;
            }
            std::ofstream results;
			std::string iteration = std::to_string(count++);
			std::string filepath;
			// default value to "test"
			if (argc > 1) {
				filepath = argv[1];
			} else {
				filepath = "test";
			}
			filepath = filepath + strDate + iteration + extension;
			std::string finalpath = folder + slash + filepath;
            if (printResults) results.open(finalpath , std::ofstream::out | std::ofstream::trunc);
			std::cout << finalpath << std::endl;

            // WRITE THE COST ON A FILE
            if (multiLevel && !stableMethod) {
                std::ofstream costFile;
                filepath = argv[1] + std::string("cost");
                filepath = filepath + iteration + extension;
                finalpath = folder + slash + filepath;
                costFile.open(finalpath, std::ofstream::out | std::ofstream::trunc);
                costFile << cost;
                costFile.close();
            }

			// WRITE RESULTS ON FILE
            if (printResults) {
                for (auto it : mcmcPath) {
                    results << it.transpose() << "\n";
                }
            }

			std::cout << "acceptance Ratio = " << accRatio << std::endl;

            if (printResults) results.close();
		}
	}

	return 0;
}