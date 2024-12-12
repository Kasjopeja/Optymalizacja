#pragma once

#include "matrix.h"

class file_reader
{
public:
	static matrix fileToMatrix(int rows, int cols, std::string filename);
};