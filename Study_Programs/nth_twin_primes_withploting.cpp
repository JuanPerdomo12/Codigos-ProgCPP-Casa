#include <iostream>
#include <fstream>
#include <cmath>

// declaration
bool is_prime(long long n);
//nth_twin_prime
int main(int argc, char **argv) {
    int n_given = std::stoll(argv[1]);
    long long counter = 0;
    double norm = 0.0;

    // open file for writing
    std::ofstream file("twin_primes.txt");
    
    for (long k = 2; ; ++k) {
        
        if (is_prime(k) == true and is_prime(k+2) == true){
            counter ++;
            norm = sqrt(k*k + (k+2)*(k+2));

            if (counter <= n_given) {
                std::cout << "nth_twin_prime(" << counter << ") is " "(" << k << ", " << k+2 << ")\n";
                file << counter << " "
                    " "<< k << "  "
                    << k+2 << "  "
                    << norm << " "
                    << std::endl;
                
            }
        }
    }
    file.close();
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