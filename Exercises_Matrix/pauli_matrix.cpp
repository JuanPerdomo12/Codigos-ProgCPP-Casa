#include "matrix_utils.h"

int main(int argc, char **argv) {
    // read size of matrix
    if (argc != 4) {
        std::cerr << "Error. Usage:\n"
                  << argv[0] << " v1 v2 v3 \n";
        return 1;
    }
    const int v1 = std::stoi(argv[1]);
    const int v2 = std::stoi(argv[2]);
    const double v3 = std::stod(argv[3]);

    // create matrix 
    std::vector<std::complex<double>> C;

    // pauli vector matrix
    char v = ' ';
    pauli_vector_matrix(v1, v2, v3, C, v);

    return 0;
}