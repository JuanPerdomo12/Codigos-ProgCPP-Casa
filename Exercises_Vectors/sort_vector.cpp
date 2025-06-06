#include <iostream>
#include <vector>
#include <algorithm>

void print(const std::vector<double>& vec) {
    for (double x : vec) {
        std::cout << x << " ";
    }
    std::cout << "\n";
}

int main() {
    std::vector<double> data = {3.2, 1.5, -4.7, 2.1, 5.3, 0.8};
    
    std::cout << "Original data: ";
    print(data);

    std::sort(data.begin(), data.end());
    /* std::sort(data.begin(), data.end(), [](double a, double b) {
    return std::abs(a) < std::abs(b);
    });
    */

    std::cout << "Sorted data: ";
    print(data);

    return 0;
}