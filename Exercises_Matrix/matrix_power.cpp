#include "matrix_utils.h"

int main(int argc, char **argv) {
    // read size of matrix
    if (argc != 4) {
        std::cerr << "Error. Usage:\n"
                  << argv[0] << " M N \n"
                  << "M : Rows\n"
                  << "N : Columns\n";
        return 1;
    }
    const int nrows = std::stoi(argv[1]);
    const int ncols = std::stoi(argv[2]);
    const int p = std::stoi(argv[3]);

    // create matrix 
    std::vector<double> M(nrows*ncols);
    std::vector<double> MM(nrows*ncols, 0.0);

    // fill matrix
    fill_matrix(M, nrows, ncols);

    // matrix to the power
    matrix_power(M, nrows, ncols, MM, p);

    // print matrix
    print_matrix(M, nrows, ncols);
    print_matrix(MM, ncols, nrows);

    return 0;
}