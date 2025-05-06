#include <iostream>
#include <fstream>
#include "derivatives.h"

double f(double x);

int main(int argc, char **argv) {
    std::cout.precision(16);
    std::cout.setf(std::ios::scientific);

    double h = 2.0e-2;

    // open file for writing
    std::ofstream file("campo.txt");
    file.precision(16);
    file.setf(std::ios::scientific);

    for(double x = 0.0; x <= 2.5; x += 2.0e-2) {
        double exact_result = -(2.0*std::sin(x) + 2.0*x*std::cos(x));
        double richardson_central_result = richardson(x, h, f, central_diff, 2);
        file << x << " " 
                << exact_result << " "
                << richardson_central_result << " "
                << std::endl;
    }
    file.close();
    return 0;
}

double f(double x)
{
   return -(2.0*x*std::sin(x));
}