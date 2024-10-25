#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
			Xopt.fit_fun(ff, ud1, ud2);
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2]{ 0,0 };
		//Tu wpisz kod funkcji
		int i = 0;
		double x1 = x0 + d;
		double fx0 = ff(x0, NAN, NAN)(0);
		double fx1 = ff(x1, NAN, NAN)(0);
		vector<double> x_vector;
		x_vector.push_back(x0);
		x_vector.push_back(x1);

		if (fx0 == fx1) {
			p[0] = x0;
			p[1] = x1;
			return p;
		}

		if (fx1 > fx0) {
			d = -d;
			x1 = x0 + d;
			x_vector[1] = x1;
			fx1 = ff(x1, NAN, NAN)(0);
			if (fx1 >= fx0) {
				p[0] = x1;
				p[1] = x0 - d;
				return p;
			}
		}

		do {
			if (i > Nmax)
				exit(-1); //W pseudokodzie bylo return error nie wiem jak to interpretowac
			i++;
			x_vector.push_back(x0 + (pow(alpha, i) * d));
		} while (ff(x_vector[i], NAN, NAN) >= ff(x_vector[i+1], NAN, NAN));

		if (d > 0) {
			p[0] = x_vector[i - 1];
			p[1] = x_vector[i + 1];
			return p;
		}

		p[0] = x_vector[i + 1];
		p[1] = x_vector[i - 1];
		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		Xopt.clear_calls;
		//Tu wpisz kod funkcji
		vector<double> ciag_fib;
		ciag_fib.push_back(1.0);
		ciag_fib.push_back(1.0);

		while (ciag_fib.back() <= (b - a) / epsilon) {
			ciag_fib.push_back(ciag_fib[ciag_fib.size() - 1] + ciag_fib[ciag_fib.size() - 2]);
		}

		int k = 1;
		while (true) {
			if (ciag_fib[k] > (b - a) / epsilon)
				break;
			k++;
		}

		double c = b - (ciag_fib[k - 1] / (ciag_fib[k] * (b - a)));
		double d = a + b - c;

		for (int i = 0; i < k - 3; i++) {
			if (ff(c, NAN, NAN)(0) < ff(d, NAN, NAN)(0)) {
				a = a; //Niepotrzebne ale spojne z pseudokodem
				b = d;
			}
			else {
				b = b; //Niepotrzebne ale spojne z pseudokodem
				a = c;
			}
			c = b - (ciag_fib[k - i - 2] / (ciag_fib[k - i - 1] * (b - a)));
			d = a + b - c;
		}

		Xopt.x = c; // Nie wiem czy uzywac tego konstruktora czy po prostu dac "= c"
		Xopt.y = ff(c, NAN, NAN)(0);
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji
		int i = 0;
		double c = (a + b) / 2;
		double d;
		double prev_d = INFINITY; //Sluzy jako d[i-1]
		while(true) {
			double l = (ff(a, NAN, NAN)(0) * (pow(b, 2) - pow(c, 2))) + (ff(b, NAN, NAN)(0) * (pow(c, 2) - pow(a, 2))) + (ff(c, NAN, NAN)(0) * (pow(a, 2) - pow(b, 2)));
			double m = (ff(a, NAN, NAN)(0) * (b - c) + ff(b, NAN, NAN)(0) * (c - a) + ff(c, NAN, NAN)(0) * (a - b));
			Xopt.f_calls += 6;
			if (m <= 0) {
				exit(-1); //W pseudokodzie bylo return error nie wiem jak to interpretowac
			}
			d = 0.5 * l / m;
			if (a < d && d < c) {
				if (ff(d, NAN, NAN)(0) < ff(c, NAN, NAN)(0)) {
					a = a;
					c = d;
					b = c;
				}
				else {
					a = d;
					c = c;
					b = b;
				}
			}
			else if (c < d && d < b) {
				if (ff(d, NAN, NAN)(0) < ff(c, NAN, NAN)(0)) {
					a = c;
					c = d;
					b = b;
				} 
				else {
					a = a;
					c = c;
					b = d;
				}
			}
			else {
				exit(-1); //W pseudokodzie bylo return error nie wiem jak to interpretowac
			}
			if (Xopt.f_calls > Nmax) {
				exit(-1); //W pseudokodzie bylo return error nie wiem jak to interpretowac
			}
			if (b - a < epsilon || abs(d - prev_d) < gamma) {
				break;
			}
			prev_d = d;
		} 

		Xopt.x = d;
		Xopt.y = ff(d, NAN, NAN)(0);

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
