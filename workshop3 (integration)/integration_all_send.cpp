#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>

//Declaration
using fptr = std::function<double(double)>;

double trapezoid(double a, double b, int N, fptr f);
double simpson(double, double b, int N, fptr f);

using algptr = std::function<double(double, double, int, fptr)>;

double richardson(double a, double b, int N, fptr f, algptr alg, int alpha);

//function
double f(double x);

int main(void){
    double exact_result = 1.56898093601156;


    // open file for writing
    std::ofstream file("integration_all.txt");
    file.precision(15);
    file.setf(std::ios::scientific);

    for (int N = 0; N <= 2000000; ){
        if (N < 10){
            N += 2;
        } else if (N < 100) {
            N += 15;
        } else if (N < 1000) {
            N += 150;
        } else if (N < 10000) {
            N += 1500;
        } else if (N < 100000) {
            N += 15000;
        } else if (N < 1000000) {
            N += 150000;
        } else if (N < 2000000) {
            N += 200000;
        } else if (N = 2000000) {
            break;
        }
    double trapezoid_result = trapezoid(-5.0, 5.0, N, f);
    double simpson_result = simpson(-5.0, 5.0, N, f);
    double richardson_trapezoid_result = richardson(-5.0, 5.0, N, f, trapezoid, 2);
    double richardson_simpson_result = richardson(-5.0, 5.0, N, f, simpson, 4);
    double error_trapezoid = std::abs(1.0 - trapezoid_result/exact_result);
    double error_simpson = std::abs(1.0 - simpson_result/exact_result);
    double error_richardson_trapezoid = std::abs(1.0 - richardson_trapezoid_result/exact_result);
    double error_richardson_simpson = std::abs(1.0 - richardson_simpson_result/exact_result);
    file << N << " "
            << trapezoid_result << " "
            << error_trapezoid << " "
            << richardson_trapezoid_result << " "
            << error_richardson_trapezoid << " "
            << simpson_result << " "
            << error_simpson << " "
            << richardson_simpson_result << " "
            << error_richardson_simpson <<
            std::endl;
    }
    file.close();
    return 0;
}

double f( double x){
    return std::sin(std::exp(x));
}

//Implementation
double trapezoid(double a, double b, int N, fptr f){
    double deltax = (b - a)/N;
    double sum = 0.0;

    for (int k = 1; k <= N-1; k++){
        double x_k = a + k * deltax;
        sum = sum + f(x_k);
    }
    double result = deltax * (0.5 * f(a) + sum + 0.5 * f(b));
    return result;
}

double simpson(double a, double b, int N, fptr f){
    if (N % 2 != 0){
        N ++;
    }

    double deltax = (b - a)/N;
    double sum1 = 0.0;
    double sum2 = 0.0;

    for (int k = 1; k <= N * 0.5; k++){
        double x_2k1 = a + (2 * k - 1) * deltax;
        sum1 = sum1 + f(x_2k1);
    }

    for (int k = 1; k <= N * 0.5 - 1; k++){
        double x_2k = a + 2 * k * deltax;
        sum2 = sum2 + f(x_2k);
    }
    double result = deltax * (1.0/3.0) * (f(a) + 4 * sum1 + 2 * sum2 + f(b));
    return result;
}

double richardson(double a, double b, int N, fptr f, algptr alg, int alpha){
    double aux = std::pow(2.0, alpha);
    double R1 = alg(a, b, N, f);
    double R2 = alg(a, b, N * 2, f);
    double result = (aux * R2 - R1) / (aux - 1.0);
    return result;
}