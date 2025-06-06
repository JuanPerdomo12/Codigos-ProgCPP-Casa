#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iomanip>

std::vector<double> read_data(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<double> data;
    double value;
    while (file >> value) {
        data.push_back(value);
    }
    return data;
}

std::vector<size_t> compute_histogram(const std::vector<double>& data, size_t bins, double min, double max) {
    std::vector<size_t> histogram(bins, 0);
    double bin_width = (max - min) / bins;

    for (double value : data) {
        if (value < min || value > max) continue;
        size_t bin = static_cast<size_t>((value - min) / bin_width);
        if (bin >= bins) bin = bins - 1;
        histogram[bin]++;
    }

    return histogram;
}

void print_results(const std::vector<size_t>& histogram, double min, double max) {
    size_t max_count = *std::max_element(histogram.begin(), histogram.end());
    double bin_width = (max - min) / histogram.size();

    for (size_t i = 0; i < histogram.size(); ++i) {
        double bin_start = min + i * bin_width;
        double bin_end = bin_start + bin_width;
        size_t bar_length = histogram[i] * 50 / max_count;

        std::cout << std::fixed << std::setprecision(2);
        std::cout << "[" << bin_start << " - " << bin_end << "): ";
        std::cout << histogram[i] << " |" << std::string(bar_length, '#') << "\n";
    }
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " filename bins min max\n";
        return 1;
    }

    std::string filename = argv[1];
    size_t bins = std::stoul(argv[2]);
    double min = std::stod(argv[3]);
    double max = std::stod(argv[4]);

    std::vector<double> data = read_data(filename);
    std::vector<size_t> histogram = compute_histogram(data, bins, min, max);

    print_results(histogram, min, max);

    return 0;
}