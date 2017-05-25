#include "Basics.hpp"

class Box {
private:
    // Coordinate of the "most bottom-left" (smaller x,y,z) vertex
    VectorXd vertCoord;

    // Length of the sides
    VectorXd r;

    // A uniform distribution density
    std::vector<std::uniform_real_distribution<double>> unifDist;

    // Dimension of the box
    unsigned int dim;

public:

    Box() {};

    Box(VectorXd& coordinates, VectorXd& radius);

    // Returns the coordinates of the box
    VectorXd getCoord(void);

    // Returns the radius
    VectorXd getRadius(void);

    // Create a Monte Carlo grid
    VectorXd randomPoint(std::default_random_engine* generator);
};