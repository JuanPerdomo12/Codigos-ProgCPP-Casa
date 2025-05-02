#include <iostream>
#include <cmath>
#include <string>

// declaration
bool is_palindrom(long long n);

int main(int argc, char **argv) {
    long long n_given = std::stoll(argv[1]);

    long long largest_palindrom = 0;

    for (long long j = n_given; j >= 0; j--) {
        for (long long k = n_given; k >= 0; k--) {
            long long product = j * k;
        if (is_palindrom(product) == true && product > largest_palindrom) {
            largest_palindrom = product;
            break;
        }
    }
}
    std::cout << "Largest palindromic product of numbers until " << n_given << " is: " << largest_palindrom << "\n";

    return 0;
}

// implementation
bool is_palindrom(long long n) {
    // supppose it is palindrom
    bool flag = true;
    
    // find palindroms
    std::string p = std::to_string(n);
    int left = 0;        
    int right = p.length() - 1;
    for (left = 0; left < right; left++, right--) {
            if (p[left] != p[right]) {
                flag = false;
                break;
            }
        } 
    return flag;
}