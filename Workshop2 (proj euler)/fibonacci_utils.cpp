#include "fibonacci_utils.h"

// implementation
long fib(long n) {
    if (n == 1) {
        return 1;

    }

    long a = 1;
    long b = 2;
    long c = 0;

    //find fibonacci numbers
    for (long k = 3; k <= n; k++) {
        c = a + b;
        a = b;
        b = c;
        
    }
    return b;
}