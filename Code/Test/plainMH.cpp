#include "mcmcTools.hpp"
#include <fstream>

int main()
{
    VectorXd firstGuess(2);
    firstGuess(0) = 1.5;
    firstGuess(1) = -1.0;

    std::vector<double> gamma = {0.01, 0.5, 2.0};
    int count = 0;

    for (auto it : gamma) {
        // Obtain results
        std::vector<VectorXd> MH, RAM;
        int nMCMC = 5000;
        double accRatioMH = 0.0, accRatioRAM = 0.0;
        MH = testMetropolis(firstGuess, nMCMC, &accRatioMH, it, false, 0.0);
        RAM = testMetropolis(firstGuess, nMCMC, &accRatioRAM, it, true, 0.4);


        // Write results on file
        std::fstream results, resultsRAM, accRatios;
        std::string resFile("resultsMH");
        std::string resRAM("resultsRAM");
        std::string resACC("resultsACC");
        std::string extension(".txt");
        std::string slash("/");
        std::string iteration = std::to_string(count++);
        std::string folder(DATA_PATH);
        std::string finalpath = folder + slash + resFile + "_" + iteration + extension;
        std::string finalpathRAM = folder + slash + resRAM + "_" + iteration + extension;
        std::string finalpathACC = folder + slash + resACC + "_" + iteration + extension;

        results.open(finalpath, std::ofstream::out | std::ofstream::trunc);
        resultsRAM.open(finalpathRAM, std::ofstream::out | std::ofstream::trunc);
        accRatios.open(finalpathACC, std::ofstream::out | std::ofstream::trunc);

        for (int i = 0; i < nMCMC; i++) {
            results << MH[i].transpose() << "\n";
            resultsRAM << RAM[i].transpose() << "\n";
        }
        accRatios << accRatioMH << "\t" << accRatioRAM << "\n";

        results.close();
        resultsRAM.close();
        accRatios.close();
    }

    return 0;
}