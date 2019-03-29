#include <Class.hpp>
#include <iostream>

int main(int argc, char* argv[])
{
    Class anObject(5.0);

    std::cout << "The number stored in my stupid class is: " << anObject.getN() << std::endl;

    std::cout << "A rather random number is: " << anObject.getM() << std::endl;

    return 0;
}
