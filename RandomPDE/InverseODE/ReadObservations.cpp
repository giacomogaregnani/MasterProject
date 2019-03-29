#include "ReadObservations.hpp"
#include <iostream>

void ReadObservations(std::vector<double>& tObs,
                      std::vector<VectorXd>& observations,
                      unsigned int nObs, unsigned int size,
                      std::string& filename)
{
    tObs.resize(nObs);
    observations.resize(nObs);


    std::fstream dataFile(filename, std::ofstream::in);
    for (unsigned int i = 0; i < nObs; i++) {
        dataFile >> tObs[i];
        observations[i].resize(size);

        for (unsigned int j = 0; j < size; j++) {
            dataFile >> observations[i](j);
        }
    }
}
