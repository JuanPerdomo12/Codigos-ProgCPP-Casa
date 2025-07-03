#include "particle.h"

#include <sstream>
#include <iomanip>

void Particle::print(void) {
    std::cout << mass << "\t" << rad << "\t"
              << R[0] << "\t" << R[1] << "\t" << R[2] << "\t"
              << V[0] << "\t" << V[1] << "\t" << V[2] << "\t"
              << F[0] << "\t" << F[1] << "\t" << F[2] << "\t";
}

std::string Particle::paraview_print(int iter){
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(5);
    
    
    oss << "x,y,z,rad,mass,vx,vy,vz,fx,fy,fz\n";
    
    oss << R[0] << "," << R[1] << "," << R[2] << ","
        << rad << "," << mass << ","
        << V[0] << "," << V[1] << "," << V[2] << ","
        << F[0] << "," << F[1] << "," << F[2];
    
    return oss.str();
}