﻿#include"opt_alg.h"

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
		double* p = new double[2] { 0, 0 };
		//Tu wpisz kod funkcji
	
		/*
		Wersja bez f_calls
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
		} while (ff(x_vector[i], NAN, NAN) >= ff(x_vector[i + 1], NAN, NAN));

		if (d > 0) {
			p[0] = x_vector[i - 1];
			p[1] = x_vector[i + 1];
			return p;
		}

		p[0] = x_vector[i + 1];
		p[1] = x_vector[i - 1];
		return p;
		*/

		//Wersja pod f_calls

		solution::clear_calls();
		solution X0(x0);
		solution X1(x0 + d);
		//solution X2;
		X0.fit_fun(ff);
		X1.fit_fun(ff);
		vector<solution> x_vector;
		x_vector.push_back(X0);
		x_vector.push_back(X1);

		if (X0.y(0) == X1.y(0)) {
			p[0] = X0.x(0);
			p[1] = X1.x(0);
			return p;
		}

		if (X1.y(0) > X0.y(0)) {
			d = -d;
			X1.x(0) = X0.x(0) + d;
			X1.fit_fun(ff);
			x_vector[1] = X1;
			if (X1.y(0) >= X0.y(0)) {
				p[0] = X1.x(0);
				p[1] = X0.x(0) - d;
				return p;
			}
		}

		int i = 0;
		do
		{
			if (X0.f_calls > Nmax)
				break;
			i++;
			x_vector.push_back(x0 + (pow(alpha, i) * d));
			x_vector[i + 1].fit_fun(ff); // Obliczenie y ostatniego elementu wektora
		} while (x_vector[i].y(0) >= x_vector[i + 1].y(0));
		if (d > 0)
		{
			p[0] = x_vector[i - 1].x(0);
			p[1] = x_vector[i + 1].x(0);
		}
		else {
			p[0] = x_vector[i + 1].x(0);
			p[1] = x_vector[i - 1].x(0);
		}
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
		Xopt.clear_calls();
		//Tu wpisz kod funkcji
		
		/*
		//Wersja bez f_calls
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

		double c = b - (ciag_fib[k - 1] / ciag_fib[k]) * (b - a);
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
			c = b - (ciag_fib[k - i - 2] / ciag_fib[k - i - 1]) * (b - a);
			d = a + b - c;
		}

		Xopt.x = c; // Nie wiem czy uzywac tego konstruktora czy po prostu dac "= c"
		Xopt.y = ff(c, NAN, NAN)(0);
		return Xopt;
		*/

		//Wersja z f_calls
		solution A(a), B(b), C, D;
		vector<double> ciag_fib;
		ciag_fib.push_back(1.0);
		ciag_fib.push_back(1.0);

		while (ciag_fib.back() <= (B.x(0) - A.x(0)) / epsilon) {
			ciag_fib.push_back(ciag_fib[ciag_fib.size() - 1] + ciag_fib[ciag_fib.size() - 2]);
		}

		int k = 1;
		while (true) {
			if (ciag_fib[k] > (B.x(0) - A.x(0)) / epsilon)
				break;
			k++;
		}

		C.x(0) = B.x(0) - (ciag_fib[k - 1] / ciag_fib[k]) * (B.x(0) - A.x(0));
		D.x(0) = A.x(0) + B.x(0) - C.x(0);
		C.fit_fun(ff);
		D.fit_fun(ff);

		for (int i = 0; i < k - 3; i++) {
			//cout << "A: " << A.x << " B: " << B.x << endl;
			if (C.y(0) < D.y(0)) {
				B.x = D.x;
			}
			else {
				A.x = C.x;
			}
			C.x(0) = B.x(0) - (ciag_fib[k - i - 2] / ciag_fib[k - i - 1]) * (B.x(0) - A.x(0));
			D.x(0) = A.x(0) + B.x(0) - C.x(0);
			C.fit_fun(ff);
			D.fit_fun(ff);
		}

		Xopt.x = C.x(0);
		Xopt.y = C.y(0);
		Xopt.flag = 1;
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
		Xopt.clear_calls();
		solution A(a), B(b), C, D, D0;
		C.x = (a + b) / 2;
		A.fit_fun(ff);
		B.fit_fun(ff);
		C.fit_fun(ff);
		double l, m;
		int i = 0;
		while (true) {
			//cout << "A: " << A.x << " B: " << B.x << endl;

			l = A.y(0) * (pow(B.x(0), 2) - pow(C.x(0), 2)) + B.y(0) * (pow(C.x(0), 2) - pow(A.x(0), 2)) + C.y(0) * (pow(A.x(0), 2) - pow(B.x(0), 2));
			m = (A.y(0) * (B.x(0) - C.x(0))) + (B.y(0) * (C.x(0) - A.x(0))) + (C.y(0) * (A.x(0) - B.x(0)));
			if (m <= 0) {
				Xopt.flag = 0;
				break;
			}

			D0.x = D.x;
			D.x = l / (2 * m);
			D.fit_fun(ff);
			if (A.x(0) < D.x(0) && D.x(0) < C.x(0))
			{
				if (D.y(0) < C.y(0))
				{
					B.x = C.x;
					C.x = D.x;
					B.fit_fun(ff);
					C.fit_fun(ff);
				}
				else
					A.x = D.x;
				A.fit_fun(ff);
			}
			else if (C.x(0) < D.x(0) && D.x(0) < B.x(0))
			{
				if (D.y(0) < C.y(0))
				{
					A.x = C.x;
					C.x = D.x;
					A.fit_fun(ff);
					C.fit_fun(ff);
				}
				else
					B.x = D.x;
				B.fit_fun(ff);
			}
			else {
				Xopt.flag = 0;
				break;
			}


			if (i > Nmax) {
				Xopt.flag = 0;
				break;
			}

			if (B.x(0) - A.x(0) < epsilon || abs(D.x(0) - D0.x(0)) <= gamma)
			{
				break;
			}
		}
		D.fit_fun(ff);
		Xopt.x = D.x(0);
		Xopt.y = D.y(0);
		Xopt.flag = 1;
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
		Xopt.clear_calls();
		solution XB, XB_temp, X(x0);
		do {
			XB.x = X.x;
			XB.fit_fun(ff);
			X = HJ_trial(ff, XB, s);
			if (X.y < XB.y) {
				do {
					XB_temp = XB;
					XB = X;
					X.x = 2 * XB.x - XB_temp.x;
					X.fit_fun(ff);
					X = HJ_trial(ff, X, s);
					if (X.f_calls > Nmax) {
						Xopt.flag = 0;
						return Xopt;
					}
				} while (X.y < XB.y);
				X = XB;
			}
			else {
				s = alpha * s;
			}
			if (X.f_calls > Nmax) {
				Xopt.flag = 0;
				return Xopt;
			}
		} while (s > epsilon);

		Xopt.x = XB.x;
		Xopt.y = XB.y;
		Xopt.flag = 1;
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
		//Tu wpisz kod 
		solution X;
		int n = get_len(XB.x);
		matrix E(n, n);
		for (int j = 0; j < n; j++) {
			E(j, j) = 1.0;
		}
		for (int j = 0; j < n; j++) {
			X.x = XB.x + (s * E[j]);
			X.fit_fun(ff);
			if (X.y < XB.y) {
				XB = X;
			}
			else {
				X.x = XB.x - (s * E[j]);
				X.fit_fun(ff);
				if (X.y < XB.y) {
					XB = X;
				}
			}
		}
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
		//Tu wpisz kod funkcji
		solution::clear_calls();
		//int i = 0;
		int n = get_len(x0);
		matrix d(n, n);
		for (int j = 0; j < n; j++) {
			d(j, j) = 1.0;
		}
		matrix lamda(n, 1), p(n, 1), s(s0);
		solution XB, X_temp;
		XB.x = x0;
		XB.fit_fun(ff);

		while (true) {

			//zacznijmy szukaæ rozwi¹zania zgodnie z zadan¹ baz¹ 
			for (int j = 0; j < n; j++) {
				X_temp.x = XB.x + (s(j) * d[j]);
				X_temp.fit_fun(ff);
				if (X_temp.y(0) < XB.y(0)) {
					XB = X_temp;
					lamda(j) += s(j);
					s(j) *= alpha;
				}
				else {
					p(j) = p(j) + 1;
					s(j) *= -beta;
				}
			}

			//czy któryœ z kroków pogorszy³ sytuacje?
			bool change = true;
			for (int j = 0; j < n; j++) {
				if (p(j) == 0 || lamda(j) == 0)
				{
					change = false;
					break;
				}

			}

			//zmiana bazy je¿eli jest to konieczne
			if (change)
			{
				matrix Q(n, n), v(n, 1);
				for (int i = 0; i < n; ++i)
					for (int j = 0; j <= i; ++j)
						Q(i, j) = lamda(i);

				Q = d * Q;
				v = Q[0] / norm(Q[0]);
				d.set_col(v, 0);
				for (int i = 1; i < n; ++i)
				{
					matrix temp(n, 1);
					for (int j = 0; j < i; ++j)
						temp = temp + (trans(Q[i]) * d[j]) * d[j];
					v = Q[i] - temp;
					d.set_col(v, i);
				}
				s = s0;
				lamda = matrix(n, 1);
				p = matrix(n, 1);
			}

			//czy wynik jest wystarczaj¹co dok³adny?
			double max_s = abs(s(0));
			for (int i = 1; i < n; ++i) {
				if (max_s < abs(s(i))) {
					max_s = abs(s(i));
				}
			}
			if (max_s < epsilon || solution::f_calls > Nmax) {
				return XB;
			}
		}
		

	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution XB;
		XB.x = x0;
		XB.fit_fun(ff);

		solution XT;
		XT = XB;

		int i = 0;

		double s = 0.5; //dlugosc boku trójkata
		double alpha = 1.0; //Wspolczynnik odbicia
		double beta = 0.5; //Wspolczynnik zwzenia
		double gamma = 2.0; //Wspolczynnik ekspansji
		double delta = 0.5; //Wspolczynnik redukcji

		do
		{
			i++;
			XT.x = XB.x;
			XT = sym_NM(ff, XB.x, s, alpha, beta, gamma, delta, epsilon, Nmax, 5.0, c);
			c *= dc;

			if (solution::f_calls > Nmax)
			{
				XT.flag = 0;
				throw std::string("Maximum amount of f_calls reached!");
			}

			if (norm(XT.x - XB.x) < epsilon)
				break;

			XB = XT;

		} while (true);

		return XT;
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
		auto max = [&](std::vector<solution> sim, int i_min) -> double
			{
				double result = 0.0;
				for (int i = 0; i < sim.size(); ++i)
				{
					double normal = norm(sim[i_min].x - sim[i].x);
					if (result < normal)
						result = normal;
				}
				return result;
			};

		int n = get_len(x0);

		//Tworzenie bazy ortogonalnej
		matrix d = matrix(n, n);
		for (int i = 0; i < n; ++i)
			d(i, i) = 1.0;

		//Tworzenie simplexu i uzupelnianie go danymi
		std::vector<solution> simplex;
		simplex.resize(n + 1);
		simplex[0].x = x0;
		simplex[0].fit_fun(ff, ud1, ud2);
		for (int i = 1; i < simplex.size(); ++i)
		{
			simplex[i].x = simplex[0].x + s * d[i - 1];
			simplex[i].fit_fun(ff, ud1, ud2);
		}

		//Indeks najmniejszej wartoœci wierzcholka simplexu
		int i_min{};
		//Indeks najwiekszej wartosci wierzcholka simplexu
		int i_max{};

		while (max(simplex, i_min) >= epsilon)
		{
			//Wyznaczanie maksymalnego i minimalnego indeksu
			i_min = 0;
			i_max = 0;
			for (int i = 1; i < simplex.size(); ++i)
			{
				if (simplex[i].y < simplex[i_min].y)
					i_min = i;
				if (simplex[i].y > simplex[i_max].y)
					i_max = i;
			}

			//Wyznaczenie œrodka ciê¿koœci
			matrix simplex_CoG{};
			for (int i = 0; i < simplex.size(); ++i)
			{
				if (i == i_max)
					continue;
				simplex_CoG = simplex_CoG + simplex[i].x;
			}
			simplex_CoG = simplex_CoG / simplex.size();

			//Obliczanie wartoœci funkcji odbitego simplexu
			solution simplex_reflected{};
			simplex_reflected.x = simplex_CoG + alpha * (simplex_CoG - simplex[i_max].x);
			simplex_reflected.fit_fun(ff, ud1, ud2);

			if (simplex_reflected.y < simplex[i_min].y)
			{
				solution simplex_expansion{};
				simplex_expansion.x = simplex_CoG + gamma * (simplex_reflected.x - simplex_CoG);
				simplex_expansion.fit_fun(ff, ud1, ud2);
				if (simplex_expansion.y < simplex_reflected.y)
					simplex[i_max] = simplex_expansion;
				else
					simplex[i_max] = simplex_reflected;
			}
			else
			{
				if (simplex[i_min].y <= simplex_reflected.y && simplex_reflected.y < simplex[i_max].y)
					simplex[i_max] = simplex_reflected;
				else
				{
					solution simplex_narrowed{};
					simplex_narrowed.x = simplex_CoG + beta * (simplex[i_max].x - simplex_CoG);
					simplex_narrowed.fit_fun(ff, ud1, ud2);
					if (simplex_narrowed.y >= simplex[i_max].y)
					{
						for (int i = 0; i < simplex.size(); ++i)
						{
							if (i == i_min)
								continue;
							simplex[i].x = delta * (simplex[i].x + simplex[i_min].x);
							simplex[i].fit_fun(ff, ud1, ud2);
						}
					}
					else
						simplex[i_max] = simplex_narrowed;
				}
			}

			if (solution::f_calls > Nmax)
			{
				simplex[i_min].flag = 0;
				throw std::string("Maximum amount of f_calls reached!");
			}

		}

		//std::cout << simplex[i_min] << "\n";

		return simplex[i_min];
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
