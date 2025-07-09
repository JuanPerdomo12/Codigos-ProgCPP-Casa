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
    if (line.empty() || line[0] == '#') continue;
    if (sscanf(line.c_str(), "%d %d", &image.width, &image.height) == 2) break;
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

void Invert_colors_filter(PGMImage& image) {
    for (int i = 0; i < image.height; ++i) {
        for (int j = 0; j < image.width; ++j) {
            image.pixels[i][j] = 255 - image.pixels[i][j];
        }
    }
}

void Blur_filter(PGMImage& image) {
    std::valarray<std::valarray<uint8_t>> temp = image.pixels;

    for (int i = 1; i < image.height - 1; ++i) {
        for (int j = 1; j < image.width - 1; ++j) {
            int sum = 0;
            for (int di = -1; di <= 1; ++di) {
                for (int dj = -1; dj <= 1; ++dj) {
                    sum += temp[i + di][j + dj];
                }
            }
            image.pixels[i][j] = sum / 9;
        }
    }
}

void Gaussian3_blur_filter(PGMImage& image) {
    const float kernel[3][3] = {
        {1.0f/16, 2.0f/16, 1.0f/16},
        {2.0f/16, 4.0f/16, 2.0f/16},
        {1.0f/16, 2.0f/16, 1.0f/16}
    };

    std::valarray<std::valarray<uint8_t>> temp = image.pixels;

    for (int i = 1; i < image.height - 1; ++i) {
        for (int j = 1; j < image.width - 1; ++j) {
            float sum = 0.0f;
            for (int di = -1; di <= 1; ++di) {
                for (int dj = -1; dj <= 1; ++dj) {
                    sum += temp[i + di][j + dj] * kernel[di + 1][dj + 1];
                }
            }
            image.pixels[i][j] = static_cast<uint8_t>(sum);
        }
    }
}

void Gaussian5_blur_filter(PGMImage& image) {
    const float kernel[5][5] = {
        {1.0f/256, 4.0f/256, 6.0f/256, 4.0f/256, 1.0f/256},
        {4.0f/256, 16.0f/256, 24.0f/256, 16.0f/256, 4.0f/256},
        {6.0f/256, 24.0f/256, 36.0f/256, 24.0f/256, 6.0f/256},
        {4.0f/256, 16.0f/256, 24.0f/256, 16.0f/256, 4.0f/256},
        {1.0f/256, 4.0f/256, 6.0f/256, 4.0f/256, 1.0f/256},
    };

    std::valarray<std::valarray<uint8_t>> temp = image.pixels;

    for (int i = 2; i < image.height - 2; ++i) {
        for (int j = 2; j < image.width - 2; ++j) {
            float sum = 0.0f;
            for (int di = -2; di <= 2; ++di) {
                for (int dj = -2; dj <= 2; ++dj) {
                    sum += temp[i + di][j + dj] * kernel[di + 2][dj + 2];
                }
            }
            image.pixels[i][j] = static_cast<uint8_t>(sum);
        }
    }
}

void Sharpen_filter(PGMImage& image) {
    const int kernel[3][3] = {
        {0, -1, 0},
        {-1, 5, -1},
        {0, -1, 0}
    };

    std::valarray<std::valarray<uint8_t>> temp = image.pixels;

    for (int i = 1; i < image.height - 1; ++i) {
        for (int j = 1; j < image.width - 1; ++j) {
            float sum = 0.0f;
            for (int di = -1; di <= 1; ++di) {
                for (int dj = -1; dj <= 1; ++dj) {
                    sum += temp[i + di][j + dj] * kernel[di + 1][dj + 1];
                }
            }
            image.pixels[i][j] = static_cast<uint8_t>(std::clamp(sum, 0.0f, 255.0f));
        }
    }
}

void Unsharp_masking(PGMImage& image) {
    const float kernel[5][5] = {
        {-1.0f/256, -4.0f/256, -6.0f/256, -4.0f/256, -1.0f/256},
        {-4.0f/256, -16.0f/256, -24.0f/256, -16.0f/256, -4.0f/256},
        {-6.0f/256, -24.0f/256, 476.0f/256, -24.0f/256, -6.0f/256},
        {-4.0f/256, -16.0f/256, -24.0f/256, -16.0f/256, -4.0f/256},
        {-1.0f/256, -4.0f/256, -6.0f/256, -4.0f/256, -1.0f/256},
    };

    std::valarray<std::valarray<uint8_t>> temp = image.pixels;

    for (int i = 2; i < image.height - 2; ++i) {
        for (int j = 2; j < image.width - 2; ++j) {
            float sum = 0.0f;
            for (int di = -2; di <= 2; ++di) {
                for (int dj = -2; dj <= 2; ++dj) {
                    sum += temp[i + di][j + dj] * kernel[di + 2][dj + 2];
                }
            }
            image.pixels[i][j] = static_cast<uint8_t>(std::clamp(sum, 0.0f, 255.0f));
        }
    }
}

void Edge_filter(PGMImage& image) {
    const int kernel[3][3] = {
        {-1, -1, -1},
        {-1, 8, -1},
        {-1, -1, -1}
    };

    std::valarray<std::valarray<uint8_t>> temp = image.pixels;

    for (int i = 1; i < image.height - 1; ++i) {
        for (int j = 1; j < image.width - 1; ++j) {
            float sum = 0.0f;
            for (int di = -1; di <= 1; ++di) {
                for (int dj = -1; dj <= 1; ++dj) {
                    sum += temp[i + di][j + dj] * kernel[di + 1][dj + 1];
                }
            }
            image.pixels[i][j] = static_cast<uint8_t>(std::clamp(sum, 0.0f, 255.0f));
        }
    }
}

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "Uso: " << argv[0] << " <archivo_entrada.pgm>" << " <nombre_filtro" << std::endl;
        return 1;
    }

    std::string inputFilename = argv[1];
    std::string filterName = argv[2];
    std::string outputFilename = "output_" + filterName + "_" + inputFilename;

    if (inputFilename.size() < 4 || inputFilename.substr(inputFilename.size() - 4) != ".pgm") {
    std::cerr << "Error: El archivo debe ser .pgm" << std::endl;
    return 1;
    }

    PGMImage image;

    if (!readPGM(inputFilename, image)) {
        return 1;
    }

    if (filterName == "invert") {
        Invert_colors_filter(image);
    } else if (filterName == "blur") {
        Blur_filter(image);
    } else if (filterName == "gaussian3") {
        Gaussian3_blur_filter(image);
    } else if (filterName == "gaussian5") {
        Gaussian5_blur_filter(image);
    } else if (filterName == "sharpen") {
        Sharpen_filter(image);
    } else if (filterName == "unsharp") {
        Unsharp_masking(image);
    } else if (filterName == "edge") {
        Edge_filter(image);
    } else {
        std::cerr << "Error: Filtro no válido. Opciones: invert, blur, gaussian3, gaussian5, sharpen, unsharp" << std::endl;
        return 1;
    }

    if (!writePGM(outputFilename, image)) {
        return 1;
    }

    std::cout << "Filtro " << filterName << " aplicado. La imagen nueva esta guardada como '" << outputFilename << "'" << std::endl;
    return 0;
}