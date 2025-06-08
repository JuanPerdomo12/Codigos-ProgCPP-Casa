#include "matrix_utils.h"

int main(int argc, char **argv) {
    // read size of matrix
    if (argc != 2) {
        std::cerr << "Error. Usage:\n"
                  << argv[0] << " epsilon \n";
        return 1;
    }
    const double epsilon = std::stod(argv[1]);

    // pauli commute
    pauli_commutation(epsilon);

    return 0;
}