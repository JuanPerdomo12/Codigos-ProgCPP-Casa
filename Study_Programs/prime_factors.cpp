#include <iostream>
#include <cmath>
#include <string>

// declaration
bool isprime(long long n);

int main(int argc, char **argv){

    long long n_given = std::stoll(argv[1]);

    std::cout << "Prime factors of " << n_given << " are: ";
    for (long long n = 2; n <= n_given; n++){
        if (n_given % n == 0 && isprime(n) == true){
            std::cout << n << " ";
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