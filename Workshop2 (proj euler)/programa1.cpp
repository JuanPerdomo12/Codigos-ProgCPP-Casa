#include <iostream>
#include <cmath>
#include <string>

#include "prime_utils.h"

int main(int argc, char **argv){

    int n_given = std::stoi(argv[1]);

    long long sum = 0;

    for (long n = 2; n <= n_given; n++){
        if (isprime(n) == true){
            sum += n; 
        }
    }
    std::cout << "Sum of prime numbers less than or equal to " << n_given << " is: "<< sum << "\n";

    return 0;
}