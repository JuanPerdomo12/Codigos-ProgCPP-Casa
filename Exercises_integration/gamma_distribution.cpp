#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <functional>

using fptr = std::function<double(double)>;

double gamma_pdf(double x, double alpha, double beta) {
    if (x <= 0) {
        return 0.0;
    }
        return std::pow(beta, alpha) / std::tgamma(alpha) * std::pow(x, alpha-1) * std::exp(-beta*x);
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

void compute_gamma_cdf(void) {
    const double alpha = 7.5;
    const double beta = 1.0;
    const int N = 1000;
    const double x_min = 0.0;
    const double x_max = 20.0;
    const int points = 100;

    auto gamma = [alpha, beta](double t) { return gamma_pdf(t, alpha, beta); };

    std::ofstream file("gamma_cdf.txt");
    file.precision(16);
    file.setf(std::ios::scientific);

    for (int i = 0; i <= points; i++) {
        double x = x_min + i*(x_max-x_min)/points;
        double simpson_cdf = simpson(0.0, x, N, gamma);
        double gauss_cdf = gauss7(0.0, x, gamma);
        file << x << " " << simpson_cdf << " " << gauss_cdf << "\n";
    }
    file.close();
}

int main(void) {
    compute_gamma_cdf();
    return 0;
}