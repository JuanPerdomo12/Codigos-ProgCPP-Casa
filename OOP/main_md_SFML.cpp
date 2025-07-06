// g++ -std=c++17 -g particle.cpp main_md_SFML.cpp -IC:/Users/jujo1/OneDrive/Documentos/CODE_LIBRARIES/SFML-3.0.0-windows-gcc-14.2.0-mingw-64-bit/SFML-3.0.0/include -LC:/Users/jujo1/OneDrive/Documentos/CODE_LIBRARIES/SFML-3.0.0-windows-gcc-14.2.0-mingw-64-bit/SFML-3.0.0/lib -lsfml-graphics-3 -lsfml-window-3 -lsfml-system-3
#include "particle.h"
#include "integrator.h"
#include "collider.h"
#include "boundary.h"
#include <vector>
#include <random>

#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <cstdint>

void initial_conditions(std::vector<Particle> & particles, int seed);

int main(int argc, char **argv) {
  const int N = 200;
  std::vector<Particle> bodies;
  bodies.resize(N); // only one particle for now

  // parameters
  std::map<std::string, double> p;
  p["T0"] = 0.0;
  p["TF"] = 50.8767;
  p["DT"] = 0.001;
  p["G"] = -9.81; //-9.81;
  p["k"] = 2600.543;

  // window
  const unsigned int WIDTH = 800, HEIGHT = 600;
  const float METERS_TO_PIXELS = 100.0f; // 200 pixels per meter
  
  // Force collider
  Collider collider(p);

  // Time initialization
  TimeIntegrator integrator(p["DT"]);

  // Boundary conditions
  Boundary bc(2.345, 0.0, 0.0, 0.0, 1.0, 0.9); // RMAX, CX, CY, CZ, EN, ET

  // initial conditions and properties
  const int seed = std::stoi(argv[1]);
  initial_conditions(bodies, seed);
  collider.computeForces(bodies); // force at t = 0
  integrator.startIntegration(bodies); // start integration algorithm
  std::cout << p["T0"] << "\t";
  bodies[0].print();
  std::cout << "\n";
  std::cout << p["T0"] << "\t";
  bodies[1].print();
  std::cout << "\n";

  // --- VISUALIZATION SETUP ---
  sf::RenderWindow window(sf::VideoMode({WIDTH, HEIGHT}), "C++ Simulation & Visualization");
  window.setFramerateLimit(240);
  sf::CircleShape particle_shape(bodies[0].rad);
  particle_shape.setFillColor(sf::Color::Cyan);
  particle_shape.setOrigin(sf::Vector2f(0.f, 0.f));

  // color setup
  std::vector<sf::Color> particle_colors(N);
  
  // rad max and min
  double min_rad = bodies[0].rad, max_rad = bodies[0].rad;
  for (const auto& p : bodies) {
    if (p.rad < min_rad) min_rad = p.rad;
    if (p.rad > max_rad) max_rad = p.rad;
  }

  // fill for size
  for (int i = 0; i < N; i++) {
    double normalized = (bodies[i].rad - min_rad) / (max_rad - min_rad);
    particle_colors[i] = sf::Color(
      255 * normalized,       // R goes up
      0,                     // G cte
      255 * (1 - normalized)  // B goes down
    );
  }

  // Time iteration
  const int niter = int((p["TF"] - p["T0"])/p["DT"]);
  for(int ii = 1; ii < niter && window.isOpen(); ++ii) {
    // check close event
    while (auto event = window.pollEvent()) {
    if (event->is<sf::Event::Closed>())
        window.close();
    }
    // md simul
    collider.computeForces(bodies);
    integrator.timeStep(bodies);
    bc.apply(bodies);

    // 2. Run the C++ visualization
    window.clear(sf::Color::Black);
    int ip = 0;
    for (const auto& p : bodies) {
      float pixel_radius = p.rad * METERS_TO_PIXELS;
      // Set the radius for this specific particle, scaled to pixels
      particle_shape.setRadius(pixel_radius);
      // Center the origin so setPosition() refers to the center of the circle
      particle_shape.setOrigin(sf::Vector2f(pixel_radius, pixel_radius));
      // color
      particle_shape.setFillColor(particle_colors[ip]);

      // Convert world coordinates to screen coordinates
      sf::Vector2f screen_pos;
      screen_pos.x = (WIDTH / 2.0f) + (p.R[0] * METERS_TO_PIXELS);
      // Flip the Y-axis: (HEIGHT - world_y_in_pixels)
      screen_pos.y = (HEIGHT / 2.0f) - (p.R[2] * METERS_TO_PIXELS);
      
      particle_shape.setPosition(screen_pos);
      window.draw(particle_shape);
      ip++;
    }
    window.display();

    // // print local info
    // double time = p["T0"] + ii*p["DT"];
    // std::cout << time << "\t";
    // bodies[0].print();
    // std::cout << "\n";
  }

  return 0;
}

void initial_conditions(std::vector<Particle> & particles, int seed)
{
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> pos_dist(-1.0, 1.0);
  std::uniform_real_distribution<double> vel_dist(-5.0, 5.0);
  std::uniform_real_distribution<double> rad_dist(0.05, 0.15);
  std::uniform_real_distribution<double> mass_dist(0.1, 0.5);

  for (auto & p : particles) {
    // positions
    p.R[0] = pos_dist(gen);
    p.R[1] = pos_dist(gen);
    p.R[2] = pos_dist(gen);

    //velocities
    p.V[0] = vel_dist(gen);
    p.V[1] = vel_dist(gen);
    p.V[2] = vel_dist(gen);

    // radius and mass
    p.rad = rad_dist(gen);
    p.mass = mass_dist(gen);
  }
}