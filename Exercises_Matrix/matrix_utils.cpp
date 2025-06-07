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
    for(int ii = 0; ii < n; ii++){
        for(int jj = 0; jj < n; jj++){
            D[ii*n+jj] = C[ii*n+jj] - I_n[ii*n+jj];         
        }
    }
    for(int ii = 0; ii < n; ii++){
        for(int jj = 0; jj < n; jj++){
            if (D[ii*n+jj] > epsilon or D[ii*n+jj] < -epsilon){
            return false;
        }
    }
    }
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
