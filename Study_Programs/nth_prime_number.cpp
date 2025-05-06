#include <iostream>
#include <cmath>

// declaration
bool is_prime(long long n);

int main(int argc, char **argv) {
    int n_given = std::stoll(argv[1]);
    long long counter = 0;
    
    for (long k = 2; ; ++k) {
        if (is_prime(k) == true){
            counter ++;
            if (counter == n_given) {
                std::cout << "nth_prime_number(" << counter << ") is " << k << "\n";
                break;
            }
        }
    }
    return 0;
}


// implementation
bool is_prime(long long n) {
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