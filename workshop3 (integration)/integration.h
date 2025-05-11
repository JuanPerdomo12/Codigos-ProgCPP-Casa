#pragma once

#include <iostream>
#include <cmath>
#include <functional>

//Declaration
using fptr = std::function<double(double)>;

double trapezoid(double a, double b, int N, fptr f);
double simpson(double, double b, int N, fptr f);

using algptr = std::function<double(double, double, int, fptr)>;

double richardson(double a, double b, int N, fptr f, algptr alg, int alpha);
