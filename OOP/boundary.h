#pragma once

#include <valarray>
#include <iostream>

class Boundary{
  double RMAX_{0.0}, EN_{0.0}, ET_{0.0}; // RMAX: radius of the sphere, EN: normal restitution coefficient, ET: tangential rest. coeff
  std::valarray<double> C_{0.0, 0.0, 0.0}; // Sphere center

  public:
    Boundary(double RMAX, double CX, double CY, double CZ, double EN, double ET) {
      RMAX_ = RMAX;
      C_[0] = CX;
      C_[1] = CY;
      C_[2] = CZ;
      EN_ = EN;
      ET_ = ET;
    }

    template <class particle_array_t>
    void apply(particle_array_t & parray) {
      // applySphericalConstraint
      for (auto & p : parray) {
        std::valarray<double> Rprime = p.R - C_;
        double rprime = std::sqrt((Rprime*Rprime).sum());
        double delta = rprime + p.rad - RMAX_;
        if (delta > 0) {
            std::valarray<double> N = Rprime/rprime;

            double Vn = (p.V * N).sum();
            std::valarray<double> V_normal = Vn * N;
            std::valarray<double> V_tangent = p.V - V_normal;

            p.R = p.R - delta*N;

            p.V = -EN_ * V_normal + ET_ * V_tangent;
        }
      }
    }
};