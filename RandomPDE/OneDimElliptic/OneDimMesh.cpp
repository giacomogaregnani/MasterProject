#include <iostream>
#include "OneDimEllipticSolver.hpp"

OneDimMesh::OneDimMesh(double xMin, double xMax, double h):
        xMin(xMin),
        xMax(xMax),
        h(h)
{
    size_t nPoints = static_cast<size_t> (round((xMax - xMin) / h)) + 1;
    points = VectorXd::LinSpaced(nPoints, xMin, xMax);
    spacings = VectorXd::Ones(nPoints - 1) * h;
}

VectorXd& OneDimMesh::getPoints()
{
    return points;
}

VectorXd& OneDimMesh::getSpacings()
{
    return spacings;
}

void OneDimMesh::setPoints(VectorXd& pointsIn)
{
    points.resize(pointsIn.size());
    points = pointsIn;
    spacings.resize(points.size() - 1);

    for (size_t i = 0; i < points.size() - 1; i++) {
        spacings(i) = points(i+1) - points(i);
    }
}
