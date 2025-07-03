#pragma once

class TimeIntegrator{
    double dt{0.0};

    public:
    TimeIntegrator(double DT) { dt = DT; }

    template <class particle_array_t>
        void startIntegration(particle_array_t & parray) {
        for (auto & p : parray) {
            p.V = p.V - p.F*dt/(2*p.mass);
        }
    }

    template <class particle_array_t>
        void timeStep(particle_array_t & parray) {
        for (auto & p : parray) {
            p.V = p.V + p.F*dt/p.mass;
            p.R = p.R + p.V*dt; 
        }    
    }
};