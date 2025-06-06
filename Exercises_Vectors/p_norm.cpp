#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <numeric>

double compute_p_norm(const std::vector<double> & v, double p) {
    auto abs_power = [p](double acc, double x) {
        return acc + std::pow(std::abs(x), p);
    };

    double sum = std::accumulate(v.begin(), v.end(), 0.0, abs_power);

    return std::pow(sum, 1.0 / p);
}

int main(int argc, char **argv) {
    std::cout.precision(16);
    std::cout.setf(std::ios::scientific);

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

    // read the p from the command line
    const double p = std::stod(argv[3]);

    // compute p norm
    std::cout << "P_norm = " << compute_p_norm(v, p) << std::endl;

    return 0;
}