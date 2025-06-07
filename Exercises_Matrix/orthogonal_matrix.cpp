#include "matrix_utils.h"

int main(int argc, char **argv) {
    // read size of matrix
    if (argc != 4) {
        std::cerr << "Error. Usage:\n"
                  << argv[0] << " M N \n"
                  << "M : Rows\n"
                  << "N : Columns\n"
                  << "Epsilon \n";
        return 1;
    }
    const int Mrows = std::stoi(argv[1]);
    const int Mcols = std::stoi(argv[2]);
    const double epsilon = std::stod(argv[3]);

    // create matrix 
    std::vector<double> M(Mrows*Mcols);

    // fill matrix
    fill_matrix(M, Mrows, Mcols);

    // print matrix
    std::cout << "Matrix A: \n";
    print_matrix(M, Mrows, Mcols);

    // check orthogonality
    orthogonal_matrix(M, Mrows, Mcols, epsilon);

    return 0;
}