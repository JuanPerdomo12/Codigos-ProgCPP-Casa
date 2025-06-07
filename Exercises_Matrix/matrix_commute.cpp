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
    const int Arows = std::stoi(argv[1]);
    const int Acols = std::stoi(argv[2]);
    const int Brows = std::stoi(argv[3]);
    const int Bcols = std::stoi(argv[4]);
    const double epsilon = std::stod(argv[5]);

    if (Acols != Brows || Arows != Bcols) {
        std::cerr << "Error. Matrices must be square.\n";
        return 1;
    }

    // create matrix 
    std::vector<double> A(Arows*Acols);

    std::vector<double> B(Brows*Bcols);

    // fill matrix
    fill_matrix(A, Arows, Acols);

    fill_matrix(B, Brows, Bcols);

    // print matrix
    std::cout << "Matrix A: \n";
    print_matrix(A, Arows, Acols);
    std::cout << "Matrix B: \n";
    print_matrix(B, Brows, Bcols);

    // check commute
    matrix_commute(A, B, Arows, Acols, Brows, Bcols, epsilon);

    return 0;
}