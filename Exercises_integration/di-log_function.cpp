#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>

using fptr = std::function<double(double)>;

double integrand(double t) {
    if (std::abs(t) < 1e-12) return 1.0;
    return -std::log(1.0 - t)/t;
}

double subs_integrand(double u){
    if (std::abs(u) < 1e-12) return 1.0;
    if (std::abs(1.0 - u) < 1e-12) return 0.0;
    return -std::log(u)/(1.0 - u);
}

double simpson(double a, double b, int N, fptr f) {
    if (N % 2 != 0) {   
        N++;
    }

    double deltax = (b - a)/N;
    double sum1 = 0.0;
    double sum2 = 0.0;

    for (int k = 1; k <= N/2; k++) {
        double x_2k1 = a + (2 * k - 1) * deltax;
        sum1 += f(x_2k1);
    }
    for (int k = 1; k <= N/2-1; k++) {
        double x_2k = a + 2 * k * deltax;
        sum2 += f(x_2k);
    }
    double result = deltax * (1.0/3.0) * (f(a) + 4 * sum1 + 2 * sum2 + f(b));
    return result;
}

double dilogarithm(double x) {
    const double error = 1e-6;
    int N = 10;
    double prev_result = 0.0;
    double result = 0.0;
    
    // Adaptive Simpson's method
    if (x < 1.0){
    do {
        prev_result = result;
        N *= 2;
        result = simpson(0.0, x, N, integrand);
    } while (std::abs(result - prev_result) > error);
} else if (x = 1.0){
    do {
        prev_result = result;
        N *= 2;
        result = simpson(0.0, x, N, subs_integrand);
    } while (std::abs(result - prev_result) > error);
}
    return result;
}


void compute_dilogarithm(void) {
    const double x_min = 0.0;
    const double x_max = 1.0;
    const int points = 100;

    std::ofstream file("dilogarithm.txt");
    file.precision(15);
    file.setf(std::ios::scientific);

    for (int i = 1; i <= points; i++) {
        double x = x_min + i*(x_max-x_min)/points;
        double li2 = dilogarithm(x);
        file << x << " " << li2 << "\n";
    }
    file.close();
}

int main(void) {
    compute_dilogarithm();
    return 0;
}