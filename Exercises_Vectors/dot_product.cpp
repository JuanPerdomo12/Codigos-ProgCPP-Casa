#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

double compute_dot_product(const std::vector<double> & v1, const std::vector<double> & v2) {
    if (v1.size() != v2.size()) {
        std::cerr << "Error: Vectors must be of the same size" << std::endl;
        return 1;
    }
    return std::inner_product(v1.begin(), v1.end(), v2.begin(), 0.0);
}

int main(int argc , char **argv)
{
    // declare the vector
    std::vector<double> v1;
    std::vector<double> v2;

    // read the size from the command line
    const int N = std::stoi(argv[1]);

    // resize the vector
    v1.resize(N, 0.0);
    v2.resize(N, 0.0);

    // fill
    int ii1 = std::stoi(argv[2]);
    auto fill1 = [&ii1](double & x){x = ii1; ii1++;};
    std::for_each(v1.begin(), v1.end(), fill1);

    int ii2 = std::stoi(argv[3]);
    auto fill2 = [&ii2](double & x){x = ii2; ii2++;};
    std::for_each(v2.begin(), v2.end(), fill2);

    // print
    auto print = [](double x){std::cout << x << " ";};

    std::for_each(v1.begin(), v1.end(), print);
    std::cout << "\n";
    std::for_each(v2.begin(), v2.end(), print);
    std::cout << "\n";

    // compute the dot product
    std::cout << "Dot product = " << compute_dot_product(v1, v2) << std::endl;

    return 0;
}
