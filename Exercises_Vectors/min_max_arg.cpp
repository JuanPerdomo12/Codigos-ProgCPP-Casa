#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

double compute_min(const std::vector<double> & v) {
    if (v.empty()) {
        std::cerr << "Error: Vector is empty" << std::endl;
        return 1;
    }
    return *std::min_element(v.begin(), v.end());
}

double compute_max(const std::vector<double> & v) {
    if (v.empty()) {
        std::cerr << "Error: Vector is empty" << std::endl;
        return 1;
    }
    return *std::max_element(v.begin(), v.end());
}

int compute_argmin(const std::vector<double> & v) {
    if (v.empty()) {
        std::cerr << "Error: Vector is empty" << std::endl;
        return -1;
    }
    return std::distance(v.begin(), std::min_element(v.begin(), v.end()));
}

int compute_argmax(const std::vector<double> & v) {
    if (v.empty()) {
        std::cerr << "Error: Vector is empty" << std::endl;
        return -1;
    }
    return std::distance(v.begin(), std::max_element(v.begin(), v.end()));
}

int main(int argc, char **argv) {
    // declare the vector
    std::vector<double> v;

    // read the size from the command line
    const int N = std::stoi(argv[1]);

    // resize the vector
    v.resize(N, 0.0);

    // fill the vector with values
    int ii = std::stoi(argv[2]);
    auto fill = [&ii](double & x) { x = ii; ii++; };
    std::for_each(v.begin(), v.end(), fill);

    // print the vector
    auto print = [](double x) { std::cout << x << " "; };
    std::for_each(v.begin(), v.end(), print);
    std::cout << "\n";

    // compute min, max, argmin, argmax
    std::cout << "Min = " << compute_min(v) << std::endl;
    std::cout << "Max = " << compute_max(v) << std::endl;
    std::cout << "Argmin = " << compute_argmin(v) << std::endl;
    std::cout << "Argmax = " << compute_argmax(v) << std::endl;

    return 0;
}