#include <iostream>
#include <fstream>
#include <cmath>

#include "integration.h"

double f(double x);

int main(void){
    //double exact_result = 1.49028031410193;
    double exact_result = 1.56898093601156;


    // open file for writing
    std::ofstream file("integration.txt");
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
    double richardson_trapezoid_result = richardson(-5.0, 5.0, N, f, trapezoid, 1);
    double richardson_simpson_result = richardson(-5.0, 5.0, N, f, simpson, 3);
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