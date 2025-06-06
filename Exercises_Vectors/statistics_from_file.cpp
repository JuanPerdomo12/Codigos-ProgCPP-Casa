#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <stdexcept>
#include <iomanip>

std::vector<double> read_num_file(const std::string & filename) {
    std::ifstream file(filename);

    std::vector<double> num;
    double val;
    while (file >> val) {
        num.push_back(val);
    }

    return num;
}

double compute_percentile(const std::vector<double> & data, double percentile) {
    double index = percentile / 100.0 * (data.size() - 1);
    size_t lower = static_cast<size_t>(std::floor(index));
    size_t upper = static_cast<size_t>(std::ceil(index));

    if (lower == upper) {
        return data[lower];
    }

    // Linear interpolation
    return data[lower] + (index - lower) * (data[upper] - data[lower]);
}

void compute_statistics(const std::vector<double> & data) {
    // Make a sorted copy for percentiles
    std::vector<double> sortedData = data;
    std::sort(sortedData.begin(), sortedData.end());

    // Basic statistics
    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    double mean = sum / data.size();
    double median = compute_percentile(sortedData, 50);

    // Min/max
    auto [minIt, maxIt] = std::minmax_element(data.begin(), data.end());

    // Output formatting
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Statistics for " << data.size() << " data points:\n";
    std::cout << "Mean:    " << mean << "\n";
    std::cout << "Median:  " << median << "\n";
    std::cout << "25th %:  " << compute_percentile(sortedData, 25) << "\n";
    std::cout << "50th %:  " << median << "\n";
    std::cout << "75th %:  " << compute_percentile(sortedData, 75) << "\n";
    std::cout << "Min:     " << *minIt << "\n";
    std::cout << "Max:     " << *maxIt << "\n";
}

int main(int argc, char** argv) {
    std::vector<double> data = read_num_file(argv[1]);
    
    compute_statistics(data);

    return 0;
}