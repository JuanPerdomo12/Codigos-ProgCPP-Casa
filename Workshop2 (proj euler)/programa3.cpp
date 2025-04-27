#include <iostream>
#include <cmath>
#include <string>

#include "fibonacci_utils.h"

int main (int argc, char **argv) {
    int n_given = std::stoi(argv[1]);

    if (n_given < 2) {
        std::cerr << "Error: Please write a number greater than or equal to 2. \n";
        return 1;
    }

    long sum = 0;

    for (long n = 1; n <= n_given; n++) {
        long term = fib(n);
        if (term > n_given) {
            break;
        }
        if (term % 2 != 0) {
            sum += term;
        }
        
    }
    std::cout << "Fibonacci series sum of odd numbers until " << n_given << " is: " << sum << "\n";
    
    return 0;
}