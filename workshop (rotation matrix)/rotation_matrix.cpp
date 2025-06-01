#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

std::vector<double> rot_x_matrix(double theta){
    std::vector<double> Rx;
    Rx.resize(9, 0.0);

    Rx[0] = 1.0;
    Rx[4] = std::cos(theta);
    Rx[5] = -std::sin(theta);
    Rx[7] = std::sin(theta);
    Rx[8] = std::cos(theta);

    return Rx;
}

std::vector<double> rot_y_matrix(double theta){
    std::vector<double> Ry;
    Ry.resize(9, 0.0);

    Ry[0] = std::cos(theta);
    Ry[2] = std::sin(theta);
    Ry[4] = 1.0;
    Ry[6] = -std::sin(theta);
    Ry[8] = std::cos(theta);

    return Ry;
}

std::vector<double> rot_z_matrix(double theta){
    std::vector<double> Rz;
    Rz.resize(9, 0.0);

    Rz[0] = std::cos(theta);
    Rz[1] = -std::sin(theta);
    Rz[3] = std::sin(theta);
    Rz[4] = std::cos(theta);
    Rz[8] = 1.0;

    return Rz;
}

void rot(const std::vector<double> & R, double & vx, double & vy, double & vz) {
    double new_x = R[0]*vx + R[1]*vy + R[2]*vz;
    double new_y = R[3]*vx + R[4]*vy + R[5]*vz;
    double new_z = R[6]*vx + R[7]*vy + R[8]*vz;
    
    vx = new_x;
    vy = new_y;
    vz = new_z;
}

int main(int argc, char **argv){
    if (argc != 7) {
        std::cerr << "Error. Usage: " << argv[0] << " Must write: vx vy vz theta_x theta_y theta_z\n";
        return 1;
    }

    double vx = std::stod(argv[1]);
    double vy = std::stod(argv[2]);
    double vz = std::stod(argv[3]);

    double theta_x = std::stod(argv[4]);
    double theta_y = std::stod(argv[5]);
    double theta_z = std::stod(argv[6]);

    std::vector<double> Rx = rot_x_matrix(theta_x);
    rot(Rx, vx, vy, vz);

    std::vector<double> Ry = rot_y_matrix(theta_y);
    rot(Ry, vx, vy, vz);

    std::vector<double> Rz = rot_z_matrix(theta_z);
    rot(Rz, vx, vy, vz);

    std::cout << vx << " " << vy << " " << vz << std::endl;

    return 0;
}
