#include "matrix_utils.h"

int main(int argc, char **argv) {
    // read size of matrix
    if (argc != 3) {
        std::cerr << "Error. Usage:\n"
                  << argv[0] << " M N \n"
                  << "M : Rows\n"
                  << "N : Columns\n";
        return 1;
    }
    const int nrows = std::stoi(argv[1]);
    const int ncols = std::stoi(argv[2]);

    // create matrix 
    std::vector<double> A(nrows*ncols);
    std::vector<double> AT(nrows*ncols, 0.0);
    std::vector<double> A_AT(nrows*nrows, 0.0);

    // fill matrix
    fill_matrix(A, nrows, ncols);

    // create transpose matrix
    transpose_matrix(A, nrows, ncols, AT);

    // multiply matrix with its transpose
    matrix_matrix_multi(A, AT, A_AT, nrows, ncols, ncols, nrows);

    // print matrix
    print_matrix(A, nrows, ncols);
    print_matrix(AT, ncols, nrows);
    print_matrix(A_AT, nrows, nrows);

    return 0;
}