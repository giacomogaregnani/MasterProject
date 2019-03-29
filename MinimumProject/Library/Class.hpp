#ifndef CLASS_HPP
#define CLASS_HPP

#include "SubLibrary/FunctionCollections.hpp"

class Class {
private:
    double mNumber;

public:
    Class() {};

    Class(double num) :
            mNumber(num)
    {}

    double getN(void) const;

    double getM(void) const;
};

#endif