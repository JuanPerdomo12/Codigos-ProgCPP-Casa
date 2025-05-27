#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <numeric>
#include <iomanip>

// linear regression
void calc_data(const std::vector<double> & xdata, const std::vector<double> & ydata, double & m, double & b, double & delta_m, double & delta_b) {
    const size_t N = xdata.size();

    double Ex = std::accumulate(xdata.begin(), xdata.end(), 0.0) / N;
    double Ey = std::accumulate(ydata.begin(), ydata.end(), 0.0) / N;
    
    double Exx = std::inner_product(xdata.begin(), xdata.end(), xdata.begin(), 0.0) / N;
    double Exy = std::inner_product(xdata.begin(), xdata.end(), ydata.begin(), 0.0) / N;
    double Eyy = std::inner_product(ydata.begin(), ydata.end(), ydata.begin(), 0.0) / N;

    double denom = Exx - Ex * Ex;
    
    m = (Exy - Ex * Ey) / denom;
    b = (Exx * Ey - Ex * Exy) / denom;

    double s_sq = (N / (N - 2.0)) * (Eyy - Ey * Ey - m * m * (Exx - Ex * Ex));
    double delta_m_sq = s_sq / (N * denom);
    double delta_b_sq = delta_m_sq * Exx;

    delta_m = std::sqrt(delta_m_sq);
    delta_b = std::sqrt(delta_b_sq);
}

int main(int argc, char **argv) {
    std::string fname = argv[1];
    std::ifstream fin(fname);

    std::vector<double> xdata;
    std::vector<double> ydata;
    double x;
    double y;
    while (fin >> x >> y) {
        xdata.push_back(x);
        ydata.push_back(y);
    }
    const double e = 1.602176634e-19;

    double m;
    double b;
    double delta_m;
    double delta_b;
    calc_data(xdata, ydata, m, b, delta_m, delta_b);

    double h_measured = m * e;
    double h_expected = 6.62607015e-34;
    double h_error = delta_m * e;

    double phi = -b;
    double phi_error = delta_b;

    std::ofstream file("data.txt");
    file.precision(16);
    file.setf(std::ios::scientific);

    file << m << " " << delta_m << " " << b << " " << delta_b << " " 
         << h_measured << " " << h_error << " " 
         << phi << " " << phi_error;
    
    file.close();

    return 0;
}