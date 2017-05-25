#include "Box.hpp"

// ======================
Box::Box(VectorXd& coordinates, VectorXd& radius)
// ======================
{
    vertCoord = coordinates;
    r = radius;
    dim = static_cast<unsigned int> (coordinates.size());
    unifDist.resize(dim);
    for (unsigned int i = 0; i < dim; i++) {
        std::uniform_real_distribution<double>::param_type distrParam(vertCoord(i), vertCoord(i) + r(i));
        unifDist[i].param(distrParam);
    }
}


// ======================
VectorXd Box::getCoord()
// ======================
{
    return vertCoord;
}

// ======================
VectorXd Box::getRadius()
// ======================
{
    return r;
}

// =====================
VectorXd Box::randomPoint(std::default_random_engine *generator)
// =====================
{
    VectorXd v(dim);
    for (unsigned int i = 0; i < dim; i++) {
        v(i) = unifDist[i](*generator);
    }
    return v;
}
