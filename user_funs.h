#pragma once

#include"ode_solver.h"


matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);

matrix ff1T(matrix, matrix = NAN, matrix = NAN);
matrix ff1R(matrix, matrix = NAN, matrix = NAN);
matrix df1(double, matrix, matrix = NAN, matrix = NAN);

matrix ff2T(matrix, matrix = NAN, matrix = NAN);
matrix df2(double, matrix, matrix = NAN, matrix = NAN);
matrix ff2R(matrix, matrix = NAN, matrix = NAN);

matrix ff3T_outside(matrix, matrix = NAN, matrix = NAN);
matrix ff3T_inside(matrix, matrix = NAN, matrix = NAN);
matrix ff3R(matrix, matrix = NAN, matrix = NAN);
matrix df3(double, matrix, matrix = NAN, matrix = NAN);

matrix ff4T(matrix, matrix = NAN, matrix = NAN);
matrix gf4T(matrix, matrix = NAN, matrix = NAN);
matrix hf4T(matrix, matrix = NAN, matrix = NAN);
matrix ff4R(matrix, matrix = NAN, matrix = NAN);
matrix gf4R(matrix, matrix = NAN, matrix = NAN);
double sigmoid(matrix theta, matrix x);

matrix ff5T(matrix, matrix = NAN, matrix = NAN);
matrix ff5R(matrix, matrix = NAN, matrix = NAN);