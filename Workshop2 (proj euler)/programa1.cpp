#include <iostream>
#include <cmath>
#include <string>

// declaration
bool isprime(long long n);

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