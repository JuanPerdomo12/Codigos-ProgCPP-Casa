#include "matrix_utils.h"

void print_matrix(const std::vector<double> & M, int nrows, int ncols){
    //std::cout.setf(std::ios::scientific);
    //std::cout.precision(15);

    for (int ii = 0; ii < nrows; ++ii) {
        for (int jj = 0; jj < ncols; ++jj) {
            std::cout << M[ii*ncols + jj] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}
void fill_matrix_random(std::vector<double> & M, const int nrows, const int ncols, const int seed){
    std::mt19937 gen(seed);
    std::uniform_real_distribution<> dis(-1, 1);
    for (int ii = 0; ii < nrows; ii++){
        for (int jj = 0; jj < ncols; jj++){
            M[ii*ncols + jj] = dis(gen);
        }
    }

}

double trace_matrix(const std::vector<double> & M, int nrows, int ncols){
    std::cout.setf(std::ios::scientific);
    std::cout.precision(15);

    double sum = 0.0;

    for (int ii = 0; ii < nrows; ++ii) {
        for (int jj = 0; jj < ncols; ++jj) {
            if (ii == jj){
                sum += M[ii*ncols + jj];
            }
        }
    }
    std::cout << "Trace = " << sum << "\n";
    return sum;
}

void fill_matrix(std::vector<double> & data, int m, int n)
{
  for (int ii = 0; ii < m; ++ii) {
    for (int jj = 0; jj < n; ++jj) {
      data[ii*n + jj] = ii*n+jj; // A_(i, j) = i*n + j = id
    }
  }
}

void matrix_matrix_multi(const std::vector<double> & M1, const std::vector<double> & M2, 
            std::vector<double> & MM, int M1rows, int M1cols, int M2rows, int M2cols){
                if (M1cols != M2rows){
                    std::cerr << "Matrix 1 must have the same columns as rows of Matrix 2\n";
                    return ;
                }

                MM.assign(M1rows*M2cols, 0.0);

                for (int i = 0; i < M1rows; i++){
                    for (int j = 0; j < M2cols; j++){
                        for (int k = 0; k < M1cols; k++){
                            MM[i*M2cols + j] += M1[i*M1cols + k] * M2[k*M2cols + j];
                        }
                    }
                }
            }

void identity_matrix(std::vector<double> & I, int nrows){
    //std::cout.setf(std::ios::scientific);
    //std::cout.precision(15);

    I.resize(nrows*nrows, 0.0);

    for (int ii = 0; ii < nrows; ++ii) {
        for (int jj = 0; jj < nrows; ++jj) {
            if (ii == jj){
                I[ii*nrows+jj]=1.0;
            }
        }
    }
}

bool check_inverse(const std::vector<double> & A, const std::vector<double> & B, const std::vector<double> & I_n, double epsilon, int n){
    std::vector<double> C;
    std::vector<double> D;

    matrix_matrix_multi(A,B,C,n,n,n,n);
    matrix_matrix_multi(B,A,D,n,n,n,n);

    std::cout << "Matrix AB: \n";
    print_matrix(C, n, n);
    std::cout << "Matrix BA: \n";
    print_matrix(D, n, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double diff1 = std::abs(C[i*n + j] - I_n[i*n + j]);
            double diff2 = std::abs(D[i*n + j] - I_n[i*n + j]);
            if (diff1 > epsilon || diff2 > epsilon) {
                std::cout << "Matrices are not inverses (failed at element [" 
                          << i << "," << j << "] with differences " 
                          << diff1 << " and " << diff2 << ")\n";
                return false;
            }
        }
    }

    std::cout << "Matrices are inverses within epsilon = " << epsilon << "\n";
    return true;
}

void hilbert_matrix(std::vector<double> & H, int nrows, int ncols){
    H.resize(nrows*ncols);

    for (int ii = 0; ii < nrows; ++ii) {
        for (int jj = 0; jj < ncols; ++jj) {
            H[ii*ncols + jj] = 1.0 / (ii + jj + 1);
        }
    }
}

void transpose_matrix(const std::vector<double> & A, int nrows, int ncols,
                      std::vector<double> & AT)
{
  for (int ii = 0; ii < nrows; ++ii) {
    for (int jj = 0; jj < ncols; ++jj) {
      AT[jj*nrows + ii] = A[ii*ncols + jj];
    }
  }
}

void matrix_power(const std::vector<double> & M, int nrows, int ncols, std::vector<double> & MM, int p){
    if (nrows != ncols){
        std::cerr << "Matrix 1 must have the same columns as rows\n";
        return ;
    }

    if (p == 0) {
        identity_matrix(MM, nrows);
        return;
    }
    if (p == 1) {
        MM = M;
        return;
    }

    std::vector<double> temp(nrows * ncols, 0.0);
    std::vector<double> result = M;


    for (int i = 1; i < p; i++){
        matrix_matrix_multi(result, M, temp, nrows, ncols, nrows, ncols);
        result = temp;
    }
    MM = result;
}

bool idempotent_matrix(const std::vector<double> & M, int nrows, int ncols, std::vector<double> & MM, int p, double epsilon){
    if (nrows != ncols){
        std::cerr << "Matrix 1 must have the same columns as rows\n";
        return 1;
    }

    if (p == 1) {
        std::cout << "All matrices are idempotent to power 1 (A¹ = A)\n";
        return true;
    }

    matrix_power(M, nrows, ncols, MM, p);

    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < ncols; ++j) {
            double diff = std::abs(MM[i*ncols + j] - M[i*ncols + j]);
            if (diff > epsilon) {
                std::cout << "Matrix is not idempotent (failed at element [" 
                          << i << "," << j << "] with difference " << diff << ")\n";
                return false;
            }
        }
    }

    std::cout << "Matrix is idempotent to power " << p 
              << " within epsilon = " << epsilon << "\n";
    return true;
}

bool matrix_commute(const std::vector<double> & A, const std::vector<double> & B,
                    int Arows, int Acols, int Brows, int Bcols, double epsilon){
    std::vector<double> C;
    std::vector<double> D;

    matrix_matrix_multi(A,B,C,Arows,Acols,Brows,Bcols);
    matrix_matrix_multi(B,A,D,Brows,Bcols,Arows,Acols);

    std::cout << "Matrix AB: \n";
    print_matrix(C, Arows, Bcols);
    std::cout << "Matrix BA: \n";
    print_matrix(D, Brows, Acols);

    for (int i = 0; i < Arows; ++i) {
        for (int j = 0; j < Bcols; ++j) {
            double diff = std::abs(C[i*Arows + j] - D[i*Brows + j]);
            if (diff > epsilon) {
                std::cout << "Matrices are not commutative (failed at element [" 
                          << i << "," << j << "] with difference " 
                          << diff << ")\n";
                return false;
            }
        }
    }

    std::cout << "Matrices are commutative within epsilon = " << epsilon << "\n";
    return true;
}

bool orthogonal_matrix(const std::vector<double> & M, int nrows, int ncols,
                       double epsilon) {
    if (nrows != ncols) {
        std::cerr << "Matrix must be square for orthogonality check.\n";
        return false;
    }

    std::vector<double> MT(nrows * ncols);
    transpose_matrix(M, nrows, ncols, MT);

    std::vector<double> I_n;
    identity_matrix(I_n, nrows);

    std::vector<double> C;
    matrix_matrix_multi(M, MT, C, nrows, ncols, nrows, ncols);

    std::cout << "Matrix AAT: \n";
    print_matrix(C, nrows, ncols);

    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < ncols; ++j) {
            double diff = std::abs(C[i*ncols + j] - I_n[i*ncols + j]);
            if (diff > epsilon) {
                std::cout << "Matrix is not orthogonal (failed at element [" 
                          << i << "," << j << "] with difference " 
                          << diff << ")\n";
                return false;
            }
        }
    }

    std::cout << "Matrix is orthogonal within epsilon = " << epsilon << "\n";
    return true;
}

void print_complex_matrix(const std::vector<std::complex<double>> & C, int nrows, int ncols) {
    for (int ii = 0; ii < nrows; ++ii) {
        for (int jj = 0; jj < ncols; ++jj) {
            std::cout << C[ii*ncols + jj] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

void conjugate_transpose(const std::vector<std::complex<double>> & C, int nrows, int ncols, std::vector<std::complex<double>> & CCT) {
    CCT.resize(ncols * nrows);

    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < ncols; ++j) {
            CCT[j*nrows + i] = std::conj(C[i*ncols + j]);
        }
    }
}

bool unitary_matrix(const std::vector<std::complex<double>> & C, int nrows, int ncols, double epsilon) {
    std::vector<std::complex<double>> CCT;
    conjugate_transpose(C, nrows, ncols, CCT);

    std::cout << "Conjugate Transpose CCT: \n";
    print_complex_matrix(CCT, nrows, ncols);

    std::vector<std::complex<double>> C_CCT;
    complex_matrix_matrix_multi(C, CCT, C_CCT, nrows, ncols, nrows, ncols);
    
    std::cout << "Matrix multi with Conjugate Transpose: \n";
    print_complex_matrix(C_CCT, nrows, ncols);

    std::vector<std::complex<double>> CI_n;
    complex_identity_matrix(CI_n, nrows);

    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < ncols; ++j) {
            double diff = std::abs(C_CCT[i*nrows + j] - CI_n[i*nrows + j]);
            if (diff > epsilon) {
                std::cout << "Matrix is not unitary (failed at [" 
                          << i << "," << j << "] with difference " << diff << ")\n";
                return false;
            }
        }
    }
    std::cout << "Matrix is unitary within epsilon = " << epsilon << "\n";
    return true;
}

bool hermitian_matrix(const std::vector<std::complex<double>> & C, int nrows, int ncols, double epsilon) {
    std::vector<std::complex<double>> CCT;
    conjugate_transpose(C, nrows, ncols, CCT);

    std::cout << "Conjugate Transpose CCT: \n";
    print_complex_matrix(CCT, nrows, ncols);

    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < ncols; ++j) {
            double diff = std::abs(C[i*nrows + j] - CCT[i*nrows + j]);
            if (diff > epsilon) {
                std::cout << "Matrix is not Hermitian (failed at [" 
                          << i << "," << j << "] with difference " << diff << ")\n";
                return false;
            }
        }
    }
    std::cout << "Matrix is hermitian within epsilon = " << epsilon << "\n";
    return true;
}

void fill_complex_matrix(std::vector<std::complex<double>> & cdata, int m, int n, const int seed) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dis(-1.0, 1.0);

    cdata.resize(m * n);
    for (int ii = 0; ii < m; ++ii) {
        for (int jj = 0; jj < n; ++jj) {
            cdata[ii * n + jj] = std::complex<double>(dis(gen), dis(gen));
        }
    }
}

void complex_matrix_matrix_multi(const std::vector<std::complex<double>> & C1, const std::vector<std::complex<double>> & C2, 
            std::vector<std::complex<double>> & C_CCT, int C1rows, int C1cols, int C2rows, int C2cols){
    if (C1cols != C2rows){
        std::cerr << "Matrix 1 must have the same columns as rows of Matrix 2\n";
        return ;
    }

    C_CCT.assign(C1rows*C2cols, 0.0);

    for (int i = 0; i < C1rows; i++){
        for (int j = 0; j < C2cols; j++){
            for (int k = 0; k < C1cols; k++){
                C_CCT[i*C2cols + j] += C1[i*C1cols + k] * C2[k*C2cols + j];
            }
        }
    }
}
void complex_identity_matrix(std::vector<std::complex<double>> & I, int nrows) {
    I.resize(nrows * nrows, std::complex<double>(0.0, 0.0));
    for (int i = 0; i < nrows; ++i) {
        I[i * nrows + i] = std::complex<double>(1.0, 0.0);
    }
}
