#include "file_reader.h"
#include <fstream>

matrix file_reader::fileToMatrix(int rows, int cols, std::string filename)
{
    matrix result = matrix(rows, cols);

    ifstream file(filename);
    if (!file.is_open()) {
        throw string("Nie mozna otworzyæ pliku: " + filename);
    }

    int i = 0;
    std::string line;
    while (std::getline(file, line))
    {
        int j = 0;
        std::istringstream lineStream(line);
        std::string value;
        while (lineStream >> value) {
            result(i, j) = std::stod(value);
            ++j;
        }
        ++i;
    }

    return result;
}