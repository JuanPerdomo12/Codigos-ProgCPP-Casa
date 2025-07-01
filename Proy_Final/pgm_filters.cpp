#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <valarray>
#include <cstdint>

struct PGMImage {
    std::string magicNum;
    int width;
    int height;
    int maxVal;
    std::valarray<std::valarray<uint8_t>> pixels;
};

bool readPGM(const std::string& fname, PGMImage& image) {
    std::ifstream file(fname, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: No se pudo abrir el archivo " << fname << std::endl;
        return false;
    }

    file >> image.magicNum;
    if (image.magicNum != "P2" && image.magicNum != "P5") {
        std::cerr << "Error: Formato PGM no soportado (" << image.magicNum << ")" << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line[0] == '#') continue;
        if (!line.empty()) break;
    }

    sscanf(line.c_str(), "%d %d", &image.width, &image.height);
    file >> image.maxVal;

    image.pixels.resize(image.height);
    for (auto& row : image.pixels) {
        row.resize(image.width);
    }

    if (image.magicNum == "P2") {
        for (int i = 0; i < image.height; ++i) {
            for (int j = 0; j < image.width; ++j) {
                int pixel;
                file >> pixel;
                image.pixels[i][j] = static_cast<uint8_t>(pixel);
            }
        }
    } else {
        file.ignore();
        for (int i = 0; i < image.height; ++i) {
            file.read(reinterpret_cast<char*>(&image.pixels[i][0]), image.width);
        }
    }

    file.close();
    return true;
}

bool writePGM(const std::string& filename, const PGMImage& image) {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: No se pudo crear el archivo " << filename << std::endl;
        return false;
    }

    file << image.magicNum << "\n";
    file << image.width << " " << image.height << "\n";
    file << image.maxVal << "\n";

    if (image.magicNum == "P2") {
        for (int i = 0; i < image.height; ++i) {
            for (int j = 0; j < image.width; ++j) {
                file << static_cast<int>(image.pixels[i][j]) << " ";
            }
            file << "\n";
        }
    } else {
        for (int i = 0; i < image.height; ++i) {
            file.write(reinterpret_cast<const char*>(&image.pixels[i][0]), image.width);
        }
    }

    file.close();
    return true;
}

int main() {
    PGMImage image;

    if (!readPGM("baboon.ascii.pgm", image)) {
        return 1;
    }

    for (int i = 0; i < image.height; ++i) {
        for (int j = 0; j < image.width; ++j) {
            image.pixels[i][j] = 255 - image.pixels[i][j];
        }
    }

    if (!writePGM("output_baboon.pgm", image)) {
        return 1;
    }

    std::cout << "Imagen nueva guardada 'output_.pgm'" << std::endl;
    return 0;
}