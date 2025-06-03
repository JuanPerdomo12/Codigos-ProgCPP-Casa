#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>
#include <vector>
#include <iomanip>

using fptr = std::function<double(double)>;

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

double gauss15(double a, double b, fptr f) {
    const std::vector<double> x = {
        -0.987992518020485, -0.937273392400706, -0.848206583410427,
        -0.724417731360170, -0.570972172608539, -0.394151347077563,
        -0.201194093997435, 0.0, 0.201194093997435,
        0.394151347077563, 0.570972172608539, 0.724417731360170,
        0.848206583410427, 0.937273392400706, 0.987992518020485
    };
    
    const std::vector<double> w = {
        0.030753241996117, 0.070366047488108, 0.107159220467172,
        0.139570677926154, 0.166269205816994, 0.186161000015562,
        0.198431485327111, 0.202578241925561, 0.198431485327111,
        0.186161000015562, 0.166269205816994, 0.139570677926154,
        0.107159220467172, 0.070366047488108, 0.030753241996117
    };

    double aux1 = (b - a)/2.0;
    double aux2 = (b + a)/2.0;
    double result = 0.0;
    
    for(int i = 0; i < 15; ++i) {
        result += w[i] * f(aux1 * x[i] + aux2);
    }
    return aux1 * result;
}

void compute_arc_length(void) {
    const int N = 100;
    const double x_min = -1.0;
    const double x_max = 1.0;

    auto arc_lenght_integrand = [](double x) {
        double dy_dx = std::sinh(x);
        return std::sqrt(1 + dy_dx * dy_dx);
    };

    std::ofstream file("arc_length.txt");
    file.precision(16);
    file.setf(std::ios::scientific);

    double simpson_arc = simpson(x_min, x_max, N, arc_lenght_integrand);
    double gauss15_arc = gauss15(x_min, x_max, arc_lenght_integrand);

    double theoric = 2 * std::sinh(1);

    file << simpson_arc << " " << gauss15_arc << " " << theoric << "\n";
    
    file.close();
}

int main(void) {
    compute_arc_length();
    return 0;
}