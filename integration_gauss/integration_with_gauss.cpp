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

double gauss2(double a, double b, fptr f);

double gauss3(double a, double b, fptr f);

double gauss5(double a, double b, fptr f);

double gauss7(double a, double b, fptr f);


//function
double f(double x);

int main(void){
    double exact_result = 2.0/3.0;


    // open file for writing
    std::ofstream file("integration_all.txt");
    file.precision(15);
    file.setf(std::ios::scientific);

    for (int N = 2; N <= 10; N++){
    double trapezoid_result = trapezoid(-2.0, 2.0, N, f);
    double simpson_result = simpson(-2.0, 2.0, N, f);
    double richardson_trapezoid_result = richardson(-2.0, 2.0, N, f, trapezoid, 2);
    double richardson_simpson_result = richardson(-2.0, 2.0, N, f, simpson, 4);
    double error_trapezoid = std::abs(1.0 - trapezoid_result/exact_result);
    double error_simpson = std::abs(1.0 - simpson_result/exact_result);
    double error_richardson_trapezoid = std::abs(1.0 - richardson_trapezoid_result/exact_result);
    double error_richardson_simpson = std::abs(1.0 - richardson_simpson_result/exact_result);
    double gauss_quadrature2 = gauss2(-2.0, 2.0, f);
    double error_gauss2 = std::abs(1.0 - gauss_quadrature2/exact_result);
    double gauss_quadrature3 = gauss3(-2.0, 2.0, f);
    double error_gauss3 = std::abs(1.0 - gauss_quadrature3/exact_result);
    double gauss_quadrature5 = gauss5(-2.0, 2.0, f);
    double error_gauss5 = std::abs(1.0 - gauss_quadrature5/exact_result);
    double gauss_quadrature7 = gauss7(-2.0, 2.0, f);
    double error_gauss7 = std::abs(1.0 - gauss_quadrature7/exact_result);
    file << N << " "
            << trapezoid_result << " "
            << error_trapezoid << " "
            << richardson_trapezoid_result << " "
            << error_richardson_trapezoid << " "
            << simpson_result << " "
            << error_simpson << " "
            << richardson_simpson_result << " "
            << error_richardson_simpson << " "
            << gauss_quadrature2 << " "
            << error_gauss2 << " "
            << gauss_quadrature3 << " "
            << error_gauss3 << " "
            << gauss_quadrature5 << " "
            << error_gauss5 << " "
            << gauss_quadrature7 << " "
            << error_gauss7 << std::endl;
    }
    file.close();
    return 0;
}

double f( double x){
    return 7*x*x*x - 8*x*x + 3*x + 3;
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

double gauss2(double a, double b, fptr f){
    //define point coordinates
    double x0 = -1.0/std::sqrt(3.0);
    double x1 = +1.0/std::sqrt(3.0);

    //define weights
    double w0 = 1.0;
    double w1 = 1.0;

    //capture integral
    double result = (b - a)/2.0 * (w0 * f(((b - a) * x0)/2.0 + (a + b)/2.0) + w1 * f(((b - a) * x1)/2.0 + (a + b)/2.0));
    return result;
}

double gauss3(double a, double b, fptr f){
    // aux
    double aux1 = (b - a)/2.0;
    double aux2 = (b + a)/2.0;

    //define point coordinates
    std::vector<double> x(3);
    x[0] = -std::sqrt(3.0/5.0);
    x[1] = 0.0;
    x[2] = +std::sqrt(3.0/5.0);

    //define weights
    std::vector<double> w(3);
    w[0] = 5.0/9.0;
    w[1] = 8.0/9.0;
    w[2] = 5.0/9.0;

    //capture integral
    double result = 0.0;
    for(int k = 0; k < 3; ++k){
        result = result + w[k] * f(aux1 * x[k] + aux2);
    }
    return aux1 * result;
}

double gauss5(double a, double b, fptr f){
    // aux
    double aux1 = (b - a)/2.0;
    double aux2 = (b + a)/2.0;

    //define point coordinates
    std::vector<double> x(5);
    x[0] = -1.0/3.0*std::sqrt(5+2*std::sqrt(10.0/7.0));
    x[1] = -1.0/3.0*std::sqrt(5-2*std::sqrt(10.0/7.0));
    x[2] = 0.0;
    x[3] = +1.0/3.0*std::sqrt(5-2*std::sqrt(10.0/7.0));
    x[4] = +1.0/3.0*std::sqrt(5+2*std::sqrt(10.0/7.0));

    //define weights
    std::vector<double> w(5);
    w[0] = (322.0-13.0*std::sqrt(70.0))/900.0;
    w[1] = (322.0+13.0*std::sqrt(70.0))/900.0;
    w[2] = 128.0/225.0;
    w[3] = (322.0+13.0*std::sqrt(70.0))/900.0;
    w[4] = (322.0-13.0*std::sqrt(70.0))/900.0;

    //capture integral
    double result = 0.0;
    for(int k = 0; k < 5; ++k){
        result = result + w[k] * f(aux1 * x[k] + aux2);
    }
    return aux1 * result;
}

double gauss7(double a, double b, fptr f){
    // aux
    double aux1 = (b - a)/2.0;
    double aux2 = (b + a)/2.0;

    //define point coordinates
    std::vector<double> x(7);
    x[0] = -0.949107912342759;
    x[1] = -0.741531185599394;
    x[2] = -0.405845151377397;
    x[3] = 0.0;
    x[4] = 0.405845151377397;
    x[5] = 0.741531185599394;
    x[6] = 0.949107912342759;

    //define weights
    std::vector<double> w(7);
    w[0] = 0.129484966168870;
    w[1] = 0.279705391489277;
    w[2] = 0.381830050505119;
    w[3] = 0.417959183673469;
    w[4] = 0.381830050505119;
    w[5] = 0.279705391489277;
    w[6] = 0.129484966168870;

    //capture integral
    double result = 0.0;
    for(int k = 0; k < 7; ++k){
        result = result + w[k] * f(aux1 * x[k] + aux2);
    }
    return aux1 * result;
}