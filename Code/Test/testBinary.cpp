#include <fstream>

using namespace std;

int main(void)
{
    double a = 0.012301232109;
    fstream test(std::string(DATA_PATH) + "/binaryFile.txt", ios::binary | ios::out | ofstream::trunc);
    test.write(reinterpret_cast<char*>(&a), sizeof(double));
    return 0;
}