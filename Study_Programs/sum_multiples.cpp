#include <iostream>
#include <cmath>
#include <string>

// declaration
long long sum_multiples(long long n, long long m1, long long m2);

int main(int argc, char **argv) {
    long long n_given = std::stoll(argv[1]);
    long long m1 = std::stoll(argv[2]);
    long long m2 = std::stoll(argv[3]);

    if (n_given < 1) {
        std::cerr << "Error: Please write a number greater than or equal to 1. \n";
        return 1;
    }

    long long sum = sum_multiples(n_given, m1, m2);
    std::cout << "Sum of multiples of " << m1 << " and " << m2 << " less than " << n_given << " is: " << sum << "\n";

    return 0;
}

// implementation
long long sum_multiples(long long n, long long m1, long long m2) {
    long long sum = 0;

    for (long long k = 1; k < n; k++) {
        if (k % m1 == 0 || k % m2 == 0) {
            sum += k;
        }
    }
    return sum;
}