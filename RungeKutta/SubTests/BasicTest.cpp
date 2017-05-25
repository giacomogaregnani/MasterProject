#include <Subdivision.hpp>

int main(int argc, char* argv[])
{
    if (argc < 6) {
        throw std::invalid_argument("Number of inputs must be 5: \n"
                                    "Number of subdivisions\n"
                                    "Time step\n"
                                    "Number of samples for mapping\n"
                                    "Number of MC samples\n"
                                    "File name for output\n"
                                    "==========================");
    }

    // Forward map specs
    forwardMapStruct FM;
    odeDef ODE;
    ODE.ode = LORENZ;
    setProblem(&ODE);
    Butcher tableau(RK4, EXPLICIT, 0);
    std::vector<double> param = ODE.refParam;
    double h = std::atof(argv[2]);
    FM.tableau = tableau;
    FM.h = h;
    FM.ODE = ODE;
    FM.paramVec = param;


    // Subdivision
    VectorXd P, radius, firstPoint;
    if (ODE.ode == LORENZ) {
        P.resize(3);
        P << -30.0, -30.0, -0.0;
        radius.resize(3);
        radius << 60.0 , 60.0 , 60.0;
        firstPoint.resize(3);
        firstPoint << 0.1, 0.1, 30.1;
    }
    if (ODE.ode == FITZNAG) {
        P.resize(2);
        P << -5.0, -5.0;
        radius.resize(2);
        radius << 10.0, 10.0;
        firstPoint.resize(2);
        firstPoint << 1.836, 0.974;
    }
    unsigned int depth = static_cast<unsigned int>(std::atoi(argv[1]));

    SubdivisionAlgorithm Algo(depth, FM, P, radius,
                              firstPoint, (unsigned int) std::stoul(argv[3]),
                              (unsigned int) std::stoul(argv[4]));
    Algo.buildBoxCollection();
    Algo.computeAttractor();
    Algo.assembleMatrix();

    std::string file = DATA_PATH + std::string(argv[5]);
    Algo.printFinestLevel(file, true);

    return 0;
}