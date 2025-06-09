#include <iostream>
#include <chrono>
#include <random>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <fstream>

void multiply(const std::vector<double> & m1, const std::vector<double> & m2, std::vector<double> & m3);

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " SEED\n";
        return 1;
    }
    const int SEED = std::atoi(argv[1]);

    std::ofstream file("matmul_time.txt");
    file.precision(16);
    file.setf(std::ios::scientific);

    for (int exp = 2; exp <= 10; ++exp) {
        int N = std::pow(2, exp);  // N = 4, 8, 16, ..., 1024

        std::mt19937 gen(SEED);
        std::uniform_real_distribution<> dist(0., 1.);

        std::vector<double> A(N*N), B(N*N), C(N*N, 0.0);
        std::generate(A.begin(), A.end(), [&]() { return dist(gen); });
        std::generate(B.begin(), B.end(), [&]() { return dist(gen); });

        // time the multiplication
        auto start = std::chrono::high_resolution_clock::now();
        multiply(A, B, C);
        auto stop = std::chrono::high_resolution_clock::now();

        // use the matrix to avoid the compiler removing it
          std::cout << C[N/2] << std::endl;

        // get the elapsed time
        std::chrono::duration<double> elapsed = stop - start;

        // print results
        std::cout << N << " " << elapsed.count() << "\n";

        file << N << " " << elapsed.count() << "\n";
    }

    file.close();

    return 0;
}

void multiply(const std::vector<double> & m1, const std::vector<double> & m2, std::vector<double> & m3)
{
    const int N = std::sqrt(m1.size());
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            for (int k = 0; k < N; k++){
                m3[i*N + j] += m1[i*N + k] * m2[k*N + j];
            }
        }
    }
}
