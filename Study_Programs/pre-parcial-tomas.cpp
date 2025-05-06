#include <iostream>
#include <fstream>
#include <cmath>

int nth_twin_primes (int N);
bool isprime(long long n);

int main (void) {
    std::ofstream file("gemelos.txt");
    for (int n=1;n<=100;++n){
        int gemelo_menor = nth_twin_primes(n);
        int gemelo_mayor = nth_twin_primes(n)+2;
        double norma = std::sqrt(gemelo_menor*gemelo_menor+gemelo_mayor*gemelo_mayor);
        file << n << " " 
                << norma << std::endl;

    }
    file.close();
    return 0;

}


int nth_twin_primes (int N){
    int encontrados = 0;
    if (N<1) std::cout << "Rango no valido";
    for (int j=3; encontrados < N; ++j){
        if (isprime(j) and isprime(j+2)){
            ++encontrados;
        }
        if (encontrados == N){
            std::cout << "(" << j << ", " << j+2 << ")";
            return j;
        }
    }
    return 0;
}

bool isprime(long long n) {
    bool flag = true; //Supone es primo
    // el codigo no tiene el caso de menor o igual a 2 pues nunca se le suministrara
    for (long long ii = 2; ii <= std::sqrt(n) + 1; ii++) {
        if (n%ii == 0) { // ii is divisor
            flag = false;
            break; // end, at least one divisor
        }
    }

    return flag;
}