#pragma once
#include <iostream>
#include <vector>
#include <cassert>
#include <random>
#include <complex>
#include <cmath>
#include <iomanip>
#include <string>
#include <sstream>

void print_matrix(const std::vector<double> & M, int nrows, int ncols);
void fill_matrix_random(std::vector<double> & M, const int nrows, const int ncols, const int seed);
double trace_matrix(const std::vector<double> & M, int nrows, int ncols);
void fill_matrix(std::vector<double> & data, int m, int n);
void matrix_matrix_multi(const std::vector<double> & M1, const std::vector<double> & M2, 
            std::vector<double> & MM, int M1rows, int M1cols, int M2rows, int M2cols);
void identity_matrix( std::vector<double> & I, int nrows);
bool check_inverse(const std::vector<double> & A, const std::vector<double> & B, const std::vector<double> &I_n, double epsilon, int n);
void hilbert_matrix(std::vector<double> & H, int nrows, int ncols);
void transpose_matrix(const std::vector<double> & A, int m, int n,
                      std::vector<double> & AT);
void matrix_power(const std::vector<double> & M, int nrows, int ncols,
                    std::vector<double> & MM, int p);
bool idempotent_matrix(const std::vector<double> & M, int nrows, int ncols,
                    std::vector<double> & MM, int p, double epsilon);
bool matrix_commute(const std::vector<double> & A, const std::vector<double> & B,
                    int Arows, int Acols, int Brows, int Bcols, double epsilon);
bool orthogonal_matrix(const std::vector<double> & M, int nrows, int ncols,
                    double epsilon);
void print_complex_matrix(const std::vector<std::complex<double>> & C, int nrows, int ncols);
void conjugate_transpose(const std::vector<std::complex<double>> & C, int nrows, int ncols,
                          std::vector<std::complex<double>> & CCT);
bool unitary_matrix(const std::vector<std::complex<double>> & C, int nrows, int ncols,
                          double epsilon);
bool hermitian_matrix(const std::vector<std::complex<double>> & C, int nrows, int ncols,
                          double epsilon);
void fill_complex_matrix(std::vector<std::complex<double>> & cdata, int m, int n, const int seed);
void complex_matrix_matrix_multi(const std::vector<std::complex<double>> & C1, const std::vector<std::complex<double>> & C2, 
            std::vector<std::complex<double>> & C_CCT, int C1rows, int C1cols, int C2rows, int C2cols);
void complex_identity_matrix(std::vector<std::complex<double>> & I, int nrows);
void pauli_vector_matrix(double v1, double v2, double v3, std::vector<std::complex<double>> & P, const char smg);
bool pauli_commutation(double epsilon);
void vandermonde_matrix(std::vector<double>& V, int nrows);
void matrix_polynomial(const std::vector<double> & A, const std::vector<double> & coeffs, std::vector<double> & result, int nrows, int ncols);

