#include "matrix_utils.h"

int main(int argc, char **argv) {
    // read size of matrix
    if (argc != 6) {
        std::cerr << "Error. Usage:\n"
                  << argv[0] << " M N \n"
                  << "M : Rows\n"
                  << "N : Columns\n"
                  << "Epsilon \n";
        return 1;
    }
    const int M1rows = std::stoi(argv[1]);
    const int M1cols = std::stoi(argv[2]);
    const int M2rows = std::stoi(argv[3]);
    const int M2cols = std::stoi(argv[4]);
    const double epsilon = std::stod(argv[5]);

    // create matrix 
    std::vector<double> A(M1rows*M1cols);

    std::vector<double> B(M2rows*M2cols);

    std::vector<double> I(M1rows*M1cols);

    //std::vector<double> A = {4, 7, 2, 6};
    //std::vector<double> B = {0.6, -0.7, -0.2, 0.4};

    // fill matrix
    fill_matrix(A, M1rows, M1cols);

    fill_matrix(B, M2rows, M2cols);

    // identity
    identity_matrix(I, M1rows);

    // print matrix
    std::cout << "Matrix A: \n";
    print_matrix(A, M1rows, M1cols);
    std::cout << "Matrix B: \n";
    print_matrix(B, M2rows, M2cols);
    std::cout << "Matrix I: \n";
    print_matrix(I, M1rows, M1cols);

    // check inverse
    check_inverse(A, B, I, epsilon, M1rows);

    return 0;
}