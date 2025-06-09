#include "matrix_utils.h"

int main(int argc, char **argv) {
    // read size of matrix
    if (argc != 2) {
        std::cerr << "Error. Usage:\n"
                  << " nrows \n";
        return 1;
    }
    const int nrows = std::stoi(argv[1]);

    // create matrix
    std::vector<double> V;

    // vandermonde matrix
    vandermonde_matrix(V, nrows);
    std::cout << "Vandermonde matrix of size " << nrows << ":\n";
    print_matrix(V, nrows, nrows);

    return 0;
}