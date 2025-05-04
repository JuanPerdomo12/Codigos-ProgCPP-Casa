#include <iostream>
#include <cmath>
#include <string>

// declaration
bool isprime(long long n);

int main(int argc, char **argv){

    long long n_given = std::stoll(argv[1]);
    long long n_given_reduced = n_given; // this is to reduce the n_given so its easier and faster for the system to find the prime factors

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

// implementation
bool isprime(long long n) {
    // suppose it is prime
    bool flag = true; 

    // some cases
    if (n < 2) {
        return false;
    }
    if (n == 2) {
        return true;
    }
    if (n%2 == 0) {
        return false;
    }

    // find divisors
    for (long long k = 3; k <= std::sqrt(n) + 1; k += 2) {
        if (n%k == 0) { // k is divisor
            flag = false;
            break; // end, at least one divisor
        }
    }

    return flag;
}