﻿/*********************************************
Kod stanowi uzupełnienie materiałów do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udostępniony na licencji CC BY-SA 3.0
Autor: dr inż. Łukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include"opt_alg.h"

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		lab1();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	//system("pause");
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;
	int Nmax = 10000;
	matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
	solution opt;
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Wahadlo
	Nmax = 1000;
	epsilon = 1e-2;
	lb = 0;
	ub = 5;
	double teta_opt = 1;
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
	ofstream Sout("symulacja_lab0.csv");
	Sout << hcat(Y[0], Y[1]);
	Sout.close();
	Y[0].~matrix();
	Y[1].~matrix();
}

void lab1()
{
	double* res = new double[2] { 0, 0 };
	double x0 = 2;
	double d = 18;
	double alpha = 1.5;
	int Nmax = 10000;
	res = expansion(ff1T, x0, d, alpha, Nmax);
	cout << res[0] << endl << res[1] << endl;


	double a = 2;
	double b = 100;
	double epsilon = 0.0001;
	res = expansion(ff1T, a, abs(a - b), alpha, Nmax);
	solution wynik;
	wynik = fib(ff1T, res[0], res[1], epsilon);
	//cout << wynik << endl;

	double gamma = 0.0001;
	wynik = lag(ff1T, res[0], res[1], epsilon, gamma, Nmax);
	cout << wynik << endl;
}

void lab2()
{

}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
