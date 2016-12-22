#include <vector>
#include <iostream>

int main ()
{
    std::vector<std::vector<int>> vectors(4);
    std::vector<int> merged = {};

    // Build vectors
    for (int j = 0; j < 4; j++) {
        for (int i = 0; i < 5; i++) {
            vectors[j].push_back(i);
        }
    }

    // Print vectors
    for (auto it : vectors) {
        for (auto itInt : it) {
            std::cout << itInt << " ";
        }
        std::cout << std::endl;
    }

    // Merge vectors


    for (size_t j = 0; j < 5; j++) {
        for (size_t i = 0; i < 4; i++) {
            merged.push_back(vectors[i][j]);
        }
    }

    // Print merged vector
    for (auto it : merged) {
        std::cout << it << " ";
    }
    std::cout << std::endl;

    return 0;
}