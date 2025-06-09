#include "matrix_utils.h"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " matrix_size coeff0 coeff1 ... coeffN\n";
        return 1;
    }
    
    const int n = std::stoi(argv[1]);
    std::vector<double> coeffs;
    
    // read coefficients from command line
    for (int i = 2; i < argc; ++i) {
        coeffs.push_back(std::stod(argv[i]));
    }
    
    // create matrix A
    std::vector<double> A;
    fill_matrix_random(A, n, n, 42);
    
    std::cout << "Matrix A (" << n << "*" << n << "):\n";
    print_matrix(A, n, n);
    
    std::cout << "Polynomial coefficients: [";
    for (double c : coeffs) std::cout << c << " ";
    std::cout << "]\n";
    
    // compute polynomial
    std::vector<double> result;
    matrix_polynomial(A, coeffs, result, n, n);
    
    std::cout << "Result of polynomial evaluation:\n";
    print_matrix(result, n, n);
    
    return 0;
}