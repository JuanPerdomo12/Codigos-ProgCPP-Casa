#include "matrix_utils.h"

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "Error. Usage:\n"
                  << argv[0] << " M N \n"
                  << "M : Rows\n"
                  << "N : Columns\n";
        return 1;
    }
    const int n = std::stoi(argv[1]);
    const double epsilon = std::stod(argv[2]);

    // real Matrix
    std::vector<double> A(n * n);
    std::vector<double> A_inv;

    // fill with random values
    fill_matrix_random(A, n, n, 42);
    
    std::cout << "Real Matrix A:\n";
    print_matrix(A, n, n);

    if (invert_matrix(A, A_inv, n, epsilon)) {
        std::cout << "Inverse of A:\n";
        print_matrix(A_inv, n, n);
    } else {
        std::cout << "Matrix A is not invertible.\n";
    }

    // complex Matrix
    std::vector<std::complex<double>> C(n * n);
    std::vector<std::complex<double>> C_inv;

    // fill with random complex numbers
    fill_complex_matrix(C, n, n, 24);
    
    std::cout << "\nComplex Matrix C:\n";
    print_complex_matrix(C, n, n);

    if (invert_complex_matrix(C, C_inv, n, epsilon)) {
        std::cout << "Inverse of C:\n";
        print_complex_matrix(C_inv, n, n);
    } else {
        std::cout << "Complex matrix C is not invertible.\n";
    }

    return 0;
}

