#include <iostream>
#include <valarray>
#include <cmath>
#include <functional>
#include <fstream>

typedef std::valarray<double> state_t;

void print(const state_t & y, double time, const char *fname)
{
    std::ofstream file(fname, std::ios::app);
        file << time << " " << y[0] << " " << y[1] << std::endl;
    file.close();
}

template <class deriv_t, class system_t, class printer_t>
void solve_euler(deriv_t fderiv, system_t & s, double tinit, double tend, double dt, printer_t writer)
{
    system_t dsdt(s.size());

    for(double t = tinit; t <= tend; t += dt) {
        fderiv(s, dsdt, t);

        s = s + dt*dsdt;

        writer(s, t, "output-euler.txt");
      }
}

template <class deriv_t, class system_t, class printer_t>
void solve_heun(deriv_t fderiv, system_t & s, double tinit, double tend, double dt, printer_t writer)
{
    system_t k1(s.size());
    system_t k2(s.size());
    system_t y1(s.size());

    for(double t = tinit; t <= tend; t += dt) {
        fderiv(s, k1, t);
        y1=s + dt*k1;
        fderiv(y1, k2 , t + dt);

        s = s + (dt/2.0)*(k1+k2);

        writer(s, t, "output-heun.txt");
      }
}


int main(int argc, char **argv)
{
    if(argc != 5) {
        std::cerr << "Uso: " << argv[0] << " dt t0 tf w" << std::endl;
        return 1;
    }

    const double dt = std::stod(argv[1]);
    const double t0 = std::stod(argv[2]);
    const double tf = std::stod(argv[3]);
    const double w = std::stod(argv[4]);

    int N = 2;
    state_t y(N);

    double y0 = 1.0; // position
    double v0 = 0.0; // velocity

    // initial conditions
    y[0] = y0;
    y[1] = v0;

    auto fderiv = [w](const state_t & y, state_t & dydt, double t) {
        dydt[0] = y[1];
        dydt[1] = -w*w*y[0];
    };
    
    // reset output files
    std::ofstream("output-euler.txt").close();
    std::ofstream("output-heun.txt").close();

    solve_euler(fderiv, y, t0, tf, dt, print);

    // reset initial conditions
    y[0] = y0; 
    y[1] = v0;

    solve_heun(fderiv, y, t0, tf, dt, print);
    return 0;
}