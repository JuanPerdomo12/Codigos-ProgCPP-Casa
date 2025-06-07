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

    if (ncols != nrows){
        std::cerr << "Matrix must be square \n";
        return 1;
    }

    // create matrix 
    std::vector<double> H(nrows*ncols);

    // fill matrix
    hilbert_matrix(H, nrows, ncols);

    // print matrix
    print_matrix(H, nrows, ncols);

    // trace of hilbert matrix
    trace_matrix(H, nrows, ncols);

    return 0;
}