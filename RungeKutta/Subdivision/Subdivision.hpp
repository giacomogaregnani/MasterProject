#include <vector>
#include <memory>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <map>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

#include "Box.hpp"

// Subdivision algorithm
class SubdivisionAlgorithm {
private:

    // The forward map
    std::shared_ptr<RungeKutta> forwardMap;
    unsigned int nTimeSteps;
    unsigned int nMap;

    // For RungeKutta class reasons, I have to store the time step too
    double h;

    // Number of levels
    size_t nLevels;

    // Outer box coordinates
    VectorXd outerBoxP;
    VectorXd outerBoxR;

    // First box
    VectorXd initBox;

    // Finer level boxes
    std::vector<bool> finerLevelBoxes;
    std::vector<unsigned long> nInDim;
    VectorXd fineRadius;

    // Visited boxes
    std::map<unsigned long, Box> visitedBoxes;

    // Density matrix
    SparseMatrix<double> P;
    unsigned int nMC;

    // Problem dimension
    unsigned int dim;

public:
    SubdivisionAlgorithm(unsigned int depth, forwardMapStruct FMParam,
                         VectorXd &initCoord, VectorXd &initRadius,
                         VectorXd &initBox, unsigned int nMap, unsigned int nMC);

    void printFinestLevel(std::string &fileName, bool density);

    void buildBoxCollection(void);

    bool findBox(VectorXd queryPoint,
                 std::vector<long int>& index,
                 unsigned long* indInGrid);

    VectorXd fromIndToCoord(std::vector<long int>& index);

    void computeAttractor();

    void assembleMatrix();
};

