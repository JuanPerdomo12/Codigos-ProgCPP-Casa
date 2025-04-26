#include <iostream>
#include <cmath>
#include <string>

#include "prime_utils.h"

int main(int argc, char **argv){

    long long n_given = std::stoll(argv[1]);
    long long n_given_reduced = n_given;

    for (long long i = 3; i <= std::sqrt(n_given_reduced) + 1; i += 2){
        if (n_given_reduced % i == 0){
            n_given_reduced /= i;
        }
    }

    for (long long n = n_given_reduced; n >= 2; n--){
        if (n_given % n == 0 && isprime(n) == true){
            std::cout << "Largest prime factor of " << n_given << " is: "<< n << "\n";
            break;
        }
    }
    
    return 0;
}