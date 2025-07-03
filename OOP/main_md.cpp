#include "particle.h"
#include "integrator.h"
#include "boundary.h"
#include "collider.h"

#include <iostream>
#include <vector>
#include <string>
#include <valarray>
#include <fstream>
#include <filesystem>

void initial_conditions(std::vector<Particle> & particles);

int main(int argc, char **argv) {
  // directory for paraview .csv files
  std::string display_folder = "DISPLAY";

  std::vector<Particle> bodies;
  bodies.resize(1); // only one particle for now

  // parameters
  std::map<std::string, double> p;
  p["T0"] = 0.0;
  p["TF"] = 50.87;
  p["DT"] = 0.01;
  p["G"] = -9.81;

  // Force collider
  Collider collider(p);

  // Time initialization
  TimeIntegrator integrator(p["DT"]);

  // Boundary conditions
  Boundary bc(2.345, 0.0, 0.0, 0.0, 0.7, 0.8); // RMAX, CX, CY, CZ, EN, ET

  // initial conditions and properties
  initial_conditions(bodies);
  collider.computeForces(bodies); // force at t = 0
  integrator.startIntegration(bodies); // start integration algorithm
  std::cout << p["T0"] << "\t";
  bodies[0].print();
  std::cout << "\n";

  // Time iteration
  const int niter = int((p["TF"] - p["T0"])/p["DT"]);

  // storage check (1GB lim)
  const size_t max_files = 1000000000 / (100 * bodies.size()); // assumes 100 bytes per file
  if (niter > max_files) {
    std::cerr << "Warning: Too many iterations (" << niter 
              << ") would create too many files.\n";
    return 1;
  }

  for(int ii = 1; ii < niter; ++ii) {
    collider.computeForces(bodies);
    integrator.timeStep(bodies);
    bc.apply(bodies);
    double time = p["T0"] + ii*p["DT"];

    std::cout << time << "\t";
    bodies[0].print();
    std::cout << "\n";

    // for paraview
    std::ostringstream fname;
  fname << display_folder << "/data_" << std::setw(4) << std::setfill('0') << ii << ".csv";
  std::ofstream pvfile(fname.str());
    pvfile << bodies[0].paraview_print(ii);
    pvfile.close();
  }

  return 0;
}

void initial_conditions(std::vector<Particle> & particles)
{
  particles[0].R[2] = 0.987;  // z is upwards, x to the right
  particles[0].V[0] = 3.126; //12.987; // z is upwards, x to the right
  particles[0].V[2] = 1.085; //4.9876; //3.987; // z is upwards, x to the right
  particles[0].rad  = 0.103;
  particles[0].mass = 0.337;
}