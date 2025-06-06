#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <random>

void shuffle(std::vector <double> & vec){
    // fill the vector (iota makes sequential numbers)
    std::iota(vec.begin(), vec.end(), 0.0); 

    // Shuffle the vector
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(vec.begin(), vec.end(), g);
}

void print(const std::vector<double>& vec) {
for (double x : vec) {
    std::cout << x << " ";
    }
    std::cout << "\n";
}
    
int main(int argc, char **argv) {
    // read variables
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " N" << std::endl;
        std::cerr << "N: size of the vector" << std::endl;
        return 1;
    }
    const int N = std::stoi(argv[1]);
    
    // set the vector
    std::vector<double> data;
    data.resize(N);
    
    // shuffle the vector
    shuffle(data);
    
    // print the vector
    print(data);
    
    return 0;
}