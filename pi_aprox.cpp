#include <iostream>
#include <cmath>

//declaracion
void pi_aprox(int n, double &pi, double &error);

int main(void){
    //notacion cientifica
    std::cout.precision(16);
    std::cout.setf(std::ios::scientific);

    double pi, error;

    for (int n = 1; n<= 20; ++n){
        pi = 0.0;
        pi_aprox(n, pi, error);
        std::cout <<"n = " << n << "\n pi = " << pi << "\n error = " << error << "\n";
    }

    return 0;
}

//Implementacion
void pi_aprox(int n, double &pi, double &error)
{
    //serie para calcular pi
    double suma = 0.0;
    for (int k = 0; k <= n; ++k){
        double serie = 1.0 / pow(16, k)*(
            4.0 / (8*k + 1) - 
            2.0 / (8*k + 4) - 
            1.0 / (8*k + 5) - 
            1.0 / (8*k + 6)
        );
        pi += serie;
    }
    //calculo de error
    error = fabs(1-pi/M_PI);
}