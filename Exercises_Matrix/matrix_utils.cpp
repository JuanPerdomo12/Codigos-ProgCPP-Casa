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
    M.resize(nrows * ncols);
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

void pauli_vector_matrix(double v1, double v2, double v3, std::vector<std::complex<double>> & P, const char smg) {
    P.resize(4);
    P[0] = std::complex<double>(v3, 0.0);
    P[1] = std::complex<double>(v1, -1.0 * v2); 
    P[2] = std::complex<double>(v1, 1.0 * v2); 
    P[3] = std::complex<double>(-v3, 0.0);

    std::cout << "Pauli vector matrix " << smg << ": \n";
    print_complex_matrix(P, 2, 2);
}

bool pauli_commutation(double epsilon) {
    std::vector<std::complex<double>> sigma_x;
    std::vector<std::complex<double>> sigma_y;
    std::vector<std::complex<double>> sigma_z;

    char x = 'x';
    char y = 'y';
    char z = 'z';

    pauli_vector_matrix(1.0, 0.0, 0.0, sigma_x, x); 
    pauli_vector_matrix(0.0, 1.0, 0.0, sigma_y, y); 
    pauli_vector_matrix(0.0, 0.0, 1.0, sigma_z, z); 

    std::vector<std::complex<double>> xy;
    std::vector<std::complex<double>> yx;
    complex_matrix_matrix_multi(sigma_x, sigma_y, xy, 2, 2, 2, 2);
    complex_matrix_matrix_multi(sigma_y, sigma_x, yx, 2, 2, 2, 2);
    
    std::vector<std::complex<double>> commu_xy(4);
    for(int i = 0; i < 4; ++i) {
        commu_xy[i] = xy[i] - yx[i];
    }
    
    std::vector<std::complex<double>> twoi_sigma_z = sigma_z;
    for(auto& val : twoi_sigma_z) {
        val *= std::complex<double>(0.0, 2.0);
    }

    bool xy_valid = true;
    for(int i = 0; i < 4; ++i) {
        if(std::abs(commu_xy[i] - twoi_sigma_z[i]) > epsilon) {
            xy_valid = false;
            break;
        }
    }

    std::vector<std::complex<double>> yz;
    std::vector<std::complex<double>> zy;
    complex_matrix_matrix_multi(sigma_y, sigma_z, yz, 2, 2, 2, 2);
    complex_matrix_matrix_multi(sigma_z, sigma_y, zy, 2, 2, 2, 2);
    
    std::vector<std::complex<double>> commu_yz(4);
    for(int i = 0; i < 4; ++i) {
        commu_yz[i] = yz[i] - zy[i];
    }
    
    std::vector<std::complex<double>> twoi_sigma_x = sigma_x;
    for(auto& val : twoi_sigma_x) {
        val *= std::complex<double>(0.0, 2.0);
    }

    bool yz_valid = true;
    for(int i = 0; i < 4; ++i) {
        if(std::abs(commu_yz[i] - twoi_sigma_x[i]) > epsilon) {
            yz_valid = false;
            break;
        }
    }

    std::vector<std::complex<double>> zx;
    std::vector<std::complex<double>> xz;
    complex_matrix_matrix_multi(sigma_z, sigma_x, zx, 2, 2, 2, 2);
    complex_matrix_matrix_multi(sigma_x, sigma_z, xz, 2, 2, 2, 2);
    
    std::vector<std::complex<double>> commu_zx(4);
    for(int i = 0; i < 4; ++i) {
        commu_zx[i] = zx[i] - xz[i];
    }
    
    std::vector<std::complex<double>> twoi_sigma_y = sigma_y;
    for(auto& val : twoi_sigma_y) {
        val *= std::complex<double>(0.0, 2.0);
    }

    bool zx_valid = true;
    for(int i = 0; i < 4; ++i) {
        if(std::abs(commu_zx[i] - twoi_sigma_y[i]) > epsilon) {
            zx_valid = false;
            break;
        }
    }

    // Print results
    std::cout << "\nCommutation Relation Results:\n";
    std::cout << "[sigma_1, sigma_2] = 2isigma_3: " << (xy_valid ? "PASSED" : "FAILED") << "\n";
    std::cout << "[sigma_2, sigma_3] = 2isigma_1: " << (yz_valid ? "PASSED" : "FAILED") << "\n";
    std::cout << "[sigma_3, sigma_1] = 2isigma_2: " << (zx_valid ? "PASSED" : "FAILED") << "\n";

    return xy_valid && yz_valid && zx_valid;
}

void vandermonde_matrix(std::vector<double>& V, int nrows) {
    V.resize(nrows * nrows);
    
    for (int i = 0; i < nrows; ++i) {
        double x = i + 1;
        for (int j = 0; j < nrows; ++j) {
            V[i * nrows + j] = pow(x, j);
        }
    }
}

void matrix_polynomial(const std::vector<double> & A, const std::vector<double> & coeffs, std::vector<double> & result, int nrows, int ncols) {
    if (nrows != ncols) {
        std::cerr << "Error: Matrix must be square for polynomial evaluation\n";
        return;
    }
    
    result.assign(nrows * ncols, 0.0);
    
    std::vector<double> temp(nrows * ncols, 0.0);
    std::vector<double> power(nrows * ncols, 0.0);
    
    for (int i = 0; i < nrows; ++i) {
        power[i * ncols + i] = 1.0;
    }
    
    for (size_t k = 0; k < coeffs.size(); ++k) {
        for (int i = 0; i < nrows * ncols; ++i) {
            result[i] += coeffs[k] * power[i];
        }
        
        if (k < coeffs.size() - 1) {
            matrix_matrix_multi(power, A, temp, nrows, ncols, nrows, ncols);
            power = temp;
        }
    }
}

bool invert_matrix(std::vector<double>& A, std::vector<double>& A_inv, int n, double epsilon) {
    A_inv.assign(n * n, 0.0);
    std::vector<double> augmented(2 * n * n);

    // Initialize augmented matrix [A | I]
    for (int i = 0; i < n; ++i){
        for (int j = 0; j < n; ++j){
            augmented[i * 2 * n + j] = A[i * n + j];           // A
            augmented[i * 2 * n + (j + n)] = (i == j) ? 1.0 : 0.0; // I
        }
    }

    // Gauss-Jordan elimination
    for (int i = 0; i < n; ++i){
        double pivot = augmented[i * 2 * n + i];
        if (std::abs(pivot) < epsilon) return false;

        // Normalize the pivot row
        for (int j = 0; j < 2 * n; ++j){
            augmented[i * 2 * n + j] /= pivot;
        }

        // Eliminate other rows
        for (int k = 0; k < n; ++k){
            if (k == i) continue;
            double factor = augmented[k * 2 * n + i];
            for (int j = 0; j < 2 * n; ++j){
                augmented[k * 2 * n + j] -= factor * augmented[i * 2 * n + j];
            }
        }
    }

    // Extract inverse matrix from augmented
    for (int i = 0; i < n; ++i){
        for (int j = 0; j < n; ++j){
            A_inv[i * n + j] = augmented[i * 2 * n + (j + n)];
        }
    }

    return true;
}

bool invert_complex_matrix(const std::vector<std::complex<double>>& C, std::vector<std::complex<double>>& C_inv, int n, double epsilon) {
    C_inv.assign(n * n, std::complex<double>(0.0, 0.0));
    std::vector<std::complex<double>> augmented(2 * n * n);

    for (int i = 0; i < n; ++i){
        for (int j = 0; j < n; ++j){
            augmented[i * 2 * n + j] = C[i * n + j];
            augmented[i * 2 * n + (j + n)] = (i == j) ? std::complex<double>(1.0, 0.0) : std::complex<double>(0.0, 0.0);
        }
    }

    for (int i = 0; i < n; ++i){
        std::complex<double> pivot = augmented[i * 2 * n + i];
        if (std::abs(pivot) < epsilon) return false;

        for (int j = 0; j < 2 * n; ++j){
            augmented[i * 2 * n + j] /= pivot;
        }

        for (int k = 0; k < n; ++k){
            if (k == i) continue;
            std::complex<double> factor = augmented[k * 2 * n + i];
            for (int j = 0; j < 2 * n; ++j){
                augmented[k * 2 * n + j] -= factor * augmented[i * 2 * n + j];
            }
        }
    }

    for (int i = 0; i < n; ++i){
        for (int j = 0; j < n; ++j){
            C_inv[i * n + j] = augmented[i * 2 * n + (j + n)];
        }
    }

    return true;
}
