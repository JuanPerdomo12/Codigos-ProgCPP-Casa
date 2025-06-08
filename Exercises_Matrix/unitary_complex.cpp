#include "matrix_utils.h"

int main(int argc, char **argv) {
    // read size of matrix
    if (argc != 5) {
        std::cerr << "Error. Usage:\n"
                  << argv[0] << " M N \n"
                  << "M : Rows\n"
                  << "N : Columns\n"
                  << "Epsilon \n";
        return 1;
    }
    const int nrows = std::stoi(argv[1]);
    const int ncols = std::stoi(argv[2]);
    const double epsilon = std::stod(argv[3]);
    const int seed = std::stoi(argv[4]);

    // create matrix 
    std::vector<std::complex<double>> C(nrows*ncols);

    /*std::vector<std::complex<double>> C = {
    {1.0, 0.0}, {0.0, 0.0},  
    {0.0, 0.0}, {1.0, 0.0}  
    };*/


    // fill matrix
    fill_complex_matrix(C, nrows, ncols, seed);

    // print matrix
    std::cout << "Matrix C: \n";
    print_complex_matrix(C, nrows, ncols);

    // check hermitian property
    unitary_matrix(C, nrows, ncols, epsilon);

    return 0;
}