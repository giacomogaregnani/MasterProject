#include "Subdivision.hpp"
#include <set>

// ======================
SubdivisionAlgorithm::SubdivisionAlgorithm(unsigned int depth, forwardMapStruct FMParam,
                                           VectorXd &initCoord, VectorXd &initRadius,
                                           VectorXd &firstCoord, unsigned int nMapIn, unsigned int nMCIn)
// ======================
{
    // Check that the points are good (if not, block everything)
    if (initCoord.size() != FMParam.ODE.size || initRadius.size() != FMParam.ODE.size) {
        throw std::invalid_argument("The outer box vertex and radius must have the same dimension as the ODE");
    }

    // Parameters of the forward map
    forwardMap = std::make_shared<RungeKutta> (FMParam.ODE, FMParam.paramVec, FMParam.tableau);
    h = FMParam.h;
    nMap = nMapIn;
    nMC = nMCIn;

    // Parameters of the outer box
    outerBoxP = initCoord;
    outerBoxR = initRadius;

    // Initial box
    initBox = firstCoord;

    // Chosen depth
    nLevels = depth;

    // Resize stuff
    dim = static_cast<unsigned int> (outerBoxP.size());
    nInDim.resize(dim);

    std::cout << "============================" << std::endl
              << "Algorithm initialized:" << std::endl
              << "Problem dimension = " << dim << std::endl
              << "N. of refinements = " << nLevels << std::endl
              << "Time integration step = " << h << std::endl
              << "============================" << std::endl;
}

// ======================
void SubdivisionAlgorithm::printFinestLevel(std::string &fileName, bool printDens)
// ======================
{
    std::ofstream output(fileName + std::string(".txt"), std::ofstream::out | std::ofstream::trunc);
    for (auto it : visitedBoxes) {
        output << it.second.getCoord().transpose() << "\t"
               << it.second.getRadius().transpose() << "\n";
    }
    if (printDens) {
        std::string fileNameMatrix = fileName + "Matrix.txt";
        saveMarket(P, fileNameMatrix);
    }
    output.close();
}

// ======================
void SubdivisionAlgorithm::buildBoxCollection()
// ======================
{
    // First level, only one box
    VectorXd newRadius, oldRadius = outerBoxR;

    // At the beginning, only one box is stored
    for (unsigned int i = 0; i < dim; i++) {
        nInDim[i] = 1;
    }
    unsigned long nBoxes = 1;

    // Find the dimensions of the finer level boxes (reduce only one dimension at a time)
    for (size_t i = 1; i < nLevels; i++) {
        size_t redDim = i % dim;
        newRadius = oldRadius;
        newRadius(redDim) /= 2;
        oldRadius = newRadius;
        nInDim[redDim] *= 2;
        nBoxes *= 2;
    }

    // Build the finer level
    finerLevelBoxes.resize(nBoxes);
    std::fill(finerLevelBoxes.begin(), finerLevelBoxes.end(), false);

    fineRadius = newRadius;
}

// ======================
bool SubdivisionAlgorithm::findBox(VectorXd queryPoint,
             std::vector<long int>& ind, unsigned long* indInGrid)
// ======================
{
    for (unsigned int i = 0; i < dim; i++) {
        ind[i] = static_cast<long int>(std::floor((queryPoint(i) - outerBoxP(i)) / fineRadius(i)));
    }

    bool isIn = true;
    for (unsigned int i = 0; i < dim; i++) {
        isIn = isIn && (ind[i] >= 0) && (ind[i] < nInDim[i]);
    }

    bool wasAlreadyVisited = false;
    long int multiplier;

    if (isIn) {
        *indInGrid = 0;
        // Find the index in the vectorized form
        for (unsigned int i = 0; i < dim; i++) {
            multiplier = 1;
            for (unsigned int j = 0; j < i; j++) {
                multiplier *= nInDim[j];
            }
            *indInGrid += ind[i] * multiplier;
        }
        wasAlreadyVisited = finerLevelBoxes[*indInGrid];
        finerLevelBoxes[*indInGrid] = true;
    } else {
        std::cout << "WARNING : YOU ARE GOING OUT OF THE BIG BOX. UNKNOWN BEHAVIOR." << std::endl;
    }

    return wasAlreadyVisited;
}

// ======================
VectorXd SubdivisionAlgorithm::fromIndToCoord(std::vector<long int> &ind)
// ======================
{
    VectorXd coordinates(dim);
    for (unsigned int i = 0; i < dim; i++) {
        coordinates(i) = ind[i] * fineRadius(i) + outerBoxP(i);
    }

    return coordinates;
}

// ======================
void SubdivisionAlgorithm::computeAttractor()
// ======================
{
    std::cout << "============================" << std::endl
              << "Building the attractor ..." << std::endl
              << "============================" << std::endl;

    // Stuff
    std::vector<long int> tmpInd(dim, 0);
    VectorXd tmpCoord(dim);
    size_t iterations = 1;
    unsigned long indInGrid;

    // Number of grid points
    VectorXd point(dim);
    std::default_random_engine generator{(unsigned int) time(NULL)};

    // Keep track of everything with counters
    unsigned long countTotBoxes = 1, countNewBoxes = 1;

    // Find initial box
    findBox(initBox, tmpInd, &indInGrid);
    tmpCoord = fromIndToCoord(tmpInd);
    visitedBoxes.insert(std::pair<unsigned long, Box> (indInGrid, Box(tmpCoord, fineRadius)));

    // Estimate the number of time steps (1.5 is a security factor)
    double maxRadius = fineRadius.maxCoeff();
    nTimeSteps = static_cast<unsigned int>(std::ceil(1.5 * maxRadius / h));
    std::cout << "Estimated number of time steps = " << nTimeSteps << std::endl;

    // Support structures
    std::vector<unsigned long> toExpand = {indInGrid};
    std::vector<unsigned long> toExpandNew = {};

    while (countNewBoxes != 0) {

        std::cout << "Iteration " << iterations++
                  << " n. boxes " << countTotBoxes
                  << " n. of new boxes " << countNewBoxes
                  << std::endl;

        countNewBoxes = 0;

        // Loop over the newly created boxes
        for (auto it : toExpand) {
            for (unsigned int j = 0; j < nMap; j++) {

                // Create a random point in the box and map it
                point = visitedBoxes[it].randomPoint(&generator);
                for (unsigned int k = 0; k < nTimeSteps; k++) {
                    point = forwardMap->oneStep(point, h);
                }

                // Search the mapped box, add it to toExpand if it was not already visited
                if (!findBox(point, tmpInd, &indInGrid)) {
                    tmpCoord = fromIndToCoord(tmpInd);
                    visitedBoxes.insert(std::pair<unsigned long, Box> (indInGrid, Box(tmpCoord, fineRadius)));
                    countNewBoxes++;
                    countTotBoxes++;
                    toExpandNew.push_back(indInGrid);
                }
            }
        }

        // Update the temporary list of boxes to expand
        toExpand = toExpandNew;
        toExpandNew.clear();
    }
}


// ======================
void SubdivisionAlgorithm::assembleMatrix()
// ======================
{
    std::cout << "============================" << std::endl
              << "Assembling the matrix ..." << std::endl
              << "============================" << std::endl;

    unsigned long N = visitedBoxes.size(), indInGrid;
    long int dist;
    P.resize((int) N, (int) N);
    std::vector<VectorXd> query;
    std::vector<long int> tmpInd(dim, 0);
    std::map<unsigned long, Box>::iterator itLoop, itMapped;
    std::map<std::pair<unsigned long, unsigned long>, double> triplets;
    std::map<std::pair<unsigned long, unsigned long>, double>::iterator itTriplets;
    std::pair<unsigned long, unsigned long> indicesInMatrix;
    VectorXd point(dim);

    std::default_random_engine generator{(unsigned int) time(NULL)};

    // Loop on the boxes
    itLoop = visitedBoxes.begin();
    unsigned long j = 0;

    while (itLoop != visitedBoxes.end()) {

        if (j % 200 == 0) {
            std::cout << "Looped on " << j << " boxes out of " << N << std::endl;
        }

        for (unsigned int k = 0; k < nMC; k++) {

            // Create a random point in the box and map it
            point = itLoop->second.randomPoint(&generator);
            for (unsigned long i = 0; i < nTimeSteps; i++) {
                point = forwardMap->oneStep(point, h);
            }

            // Find the box in which the mapped point is and its index in the vector
            findBox(point, tmpInd, &indInGrid);
            itMapped = visitedBoxes.find(indInGrid);

            // Security if. It shouldn't happen (will produce a small bias in the density)
            if (itMapped != visitedBoxes.end()) {

                // Find the "index" of the mapped box in the map
                dist = std::distance(visitedBoxes.begin(), itMapped);

                // Update the value of the matrix or create a new triplet
                indicesInMatrix = std::pair<unsigned long, unsigned long>(dist, j);
                itTriplets = triplets.find(indicesInMatrix);
                if (itTriplets != triplets.end()) {
                    itTriplets->second += 1.0;
                } else {
                    triplets.insert(std::pair<std::pair<unsigned long, unsigned long>,
                            double>(indicesInMatrix, 1.0));
                }
            }
        }
        ++j;
        ++itLoop;
    }

    // Convert the triplets into Eigen triplets for convenience of data storage
    std::vector<Triplet<double>> EigenTriplets;
    for (auto it : triplets) {
        EigenTriplets.push_back(Triplet<double>((int) it.first.first, (int) it.first.second, it.second));
    }
    P.setFromTriplets(EigenTriplets.begin(), EigenTriplets.end());

    // Normalize
    P /= nMC;
}