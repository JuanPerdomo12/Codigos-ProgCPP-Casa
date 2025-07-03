#pragma once

#include <iostream>
#include <valarray>
#include <string>

struct Particle {
    std::valarray<double> R{0.0, 0.0, 0.0};
    std::valarray<double> V{0.0, 0.0, 0.0};
    std::valarray<double> F{0.0, 0.0, 0.0};
    double mass{0}, rad{0};
    void print(void);
    std::string paraview_print(int iter);
};