#include <iostream>
#include <algorithm>

long long gcd(long long a, long long b) {
    while (b != 0) {
        long long temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

long long find_LCM(long n) {
    long long lcm = 1;
    for (long i = 2; i <= n; ++i) {
        lcm = (lcm * i) / gcd(lcm, i);
    }
    return lcm;
}

int main(int argc, char **argv) {
    long long n_given = std::stoll(argv[1]);

    long long result = find_LCM(n_given);

    std::cout << "The least common multiple of numbers between 1 and " << n_given << " is " << result << "\n";

    return 0;
}