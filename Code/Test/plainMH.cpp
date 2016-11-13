#include "mcmcTools.hpp"
#include <fstream>

int main()
{
    VectorXd firstGuess(2);
    firstGuess(0) = 1.5;
    firstGuess(1) = -1.0;

    std::vector<VectorXd> MH, RAM;
    int nMCMC = 5000;
    double accRatioMH, accRatioRAM;
    double gamma = 2.0;
    MH = testMetropolis(firstGuess, nMCMC, &accRatioMH, gamma, false, 0.0);
    RAM = testMetropolis(firstGuess, nMCMC, &accRatioRAM, gamma, true, 0.5);

    std::fstream results, resultsRAM;
    std::string resFile("resultsMH");
    std::string resRAM("resultsRAM");
    std::string extension(".txt");
    std::string slash("/");
    std::string folder(DATA_PATH);
    std::string finalpath = folder + slash + resFile + extension;
    std::string finalpathRAM = folder + slash + resRAM + extension;

    results.open(finalpath , std::ofstream::out | std::ofstream::trunc);
    resultsRAM.open(finalpathRAM , std::ofstream::out | std::ofstream::trunc);

    for (int i = 0; i < nMCMC; i++) {
        results << MH[i].transpose() << "\n";
        resultsRAM << RAM[i].transpose() << "\n";
    }
    results.close();

    return 0;
}