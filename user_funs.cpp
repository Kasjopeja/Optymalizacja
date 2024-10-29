#include"user_funs.h"

matrix ff0T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
	int n = get_len(Y[0]);
	double teta_max = Y[1](0, 0);
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);
	double m = 1, l = 0.5, b = 0.5, g = 9.81;
	double I = m*pow(l, 2);
	dY(0) = Y(1);
	dY(1) = ((t <= ud2(1))*ud2(0) - m*g*l*sin(Y(0)) - b*Y(1)) / I;
	return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2) {
	return -cos(0.1 * x(0)) * exp(-pow((0.1 * x(0) - 2 * 3.14), 2)) + 0.002 * pow((0.1 * x(0)), 2);
}

matrix df1(double t, matrix Y, matrix ud1, matrix ud2) {
    matrix dY(Y);  // dY to macierz przechowująca pochodne VA', VB' i TB'

    // Parametry stałe
    double a = 0.98;            // wsp. lepkości
    double b = 0.63;            // wsp. zwężenia strumienia
    double g = 9.81;            // przysp. ziemskie
    double PA = 0.5;            // pole pow. zbiornika A
    double PB = 1;              // pole pow. zbiornika B
    double DB = 0.00365665;     // wielkość otworu w zbiorniku B
    double Fin = 0.01;          // ilość wody wpływającej do B (z rury)
    double Tin = 20;            // temp. wody wpływającej
    double TA = 90.0;           // temp. wody w zbiorniku A
    double DA = ud1(0);         // wielkość otworu w zbiorniku A (parametr z ud2)

    // Stan w danym momencie, podany przez wektor Y
    double VA = Y(0);           // objętość w zbiorniku A
    double VB = Y(1);           // objętość w zbiorniku B
    double TB = Y(2);           // temperatura w zbiorniku B

    // Obliczenie przepływu wody z A i B
    double FAout = VA > 0 ? a * b * DA * sqrt(2 * g * VA / PA) : 0;
    double FBout = VB > 0 ? a * b * DB * sqrt(2 * g * VB / PB) : 0;

    // Obliczenie pochodnych VA', VB' i TB'
    dY(0) = -FAout;                                 // VA' (zmiana objętości w A)
    dY(1) = FAout + Fin - FBout;                    // VB' (zmiana objętości w B)
    dY(2) = FAout / VB * (TA - TB) + Fin / VB * (Tin - TB); // TB' (zmiana temperatury w B)

    return dY;
}

matrix ff2T(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	matrix Y0 = matrix(3, new double[3] {5, 1, 20});
	matrix* Y = solve_ode(df1, 0, 1, 2000, Y0, x);
	int n = get_len(Y[0]);
	double max = Y[1](0, 2);
	for (int i = 0; i < n; i++) {
		if (max < Y[1](i, 2)) {
			max = Y[1](i, 2);
		}
	}
	y = abs(max - 50);
	return y[0];
}


matrix ff3T(matrix x, matrix ud1, matrix ud2) {
	double x1 = x(0);
	double x2 = x(1);

	double result = pow(x1, 2) + pow(x2, 2) - cos(2.5 * 3.14 * x1) - cos(2.5 * 3.14 * x2) + 2;

	return matrix(1, 1, result);
}
