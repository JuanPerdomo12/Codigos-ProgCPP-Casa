#include <iostream>
#include <fstream>
#include <string>
#include <vector>

void readPGM(const std::string & fname) {
    std::ifstream file(fname);
    if (!file.is_open()) {
        std::cerr << "Error al abrir el archivo." << std::endl;
        return;
    }

    std::string line;
    int width;
    int height;
    int maxVal;

    // leer encabezado
    std::getline(file, line); // P2
    if (line != "P2") {
        std::cerr << "No es un PGM válido (P2)." << std::endl;
        return;
    }

    // Saltar comentarios
    while (std::getline(file, line)) {
        if (line[0] != '#') break;
    }

    // leer ancho, alto y max
    sscanf(line.c_str(), "%d %d", &width, &height);
    file >> maxVal;

    // leer pixeles
    std::vector<std::vector<int>> pixels(height, std::vector<int>(width));
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            file >> pixels[i][j];
        }
    }

    // mostrar en consola (simulación ASCII art)
    std::cout << "Imagen PGM (" << width << "x" << height << "):\n";
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            char shade = ' ';
            if (pixels[i][j] > 200) shade = '@';
            else if (pixels[i][j] > 150) shade = '#';
            else if (pixels[i][j] > 100) shade = '*';
            else if (pixels[i][j] > 50) shade = '.';
            std::cout << shade << " ";
        }
        std::cout << "\n";
    }
}

int main() {
    readPGM("apollonian_gasket.ascii.pgm");
    return 0;
}