/*********************************************
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
	std::ofstream Sout("symulacja_lab1.csv");
	
	//zadanie teoretyczne

	double* res = new double[2] { 0, 0 };
	double x0 = -100, d = 5, alpha = 2.1;
	int Nmax = 10000;

	//double a = 50, b = 70;
	double epsilon = 0.0001;
	double gamma = 0.000001;
	solution wynik;

	for (int i = 0; i < 100; i++)
	{
		res = expansion(ff1T, x0, d, alpha, Nmax);
		cout << res[0] << endl << res[1] << endl << solution::f_calls << endl << endl;
		//Sout << "x" << res[0] << ";" << "x" << res[1] << ";" << "x" << solution::f_calls << "\n";

		wynik = fib(ff1T, res[0], res[1], epsilon);
		//Sout << "x" << wynik.x << "x" << wynik.y << "x" << wynik.f_calls << "\n";
		cout << wynik << endl;

		wynik = lag(ff1T, res[0], res[1], epsilon, gamma, Nmax);
		//Sout << "x" << wynik.x << "x" << wynik.y << "x" << wynik.f_calls << "\n";
		cout << wynik << endl;

		x0 = x0 + 2;
	}



	//zadanie praktyczne

	//double* res = new double[2] { 0, 0 };
	//double da = 0.005;
	//double delta_da = 0.002;
	//double tmax = 0;
	//double alpha = 1.5;
	//double epsilon = 0.0001;
	//double gamma = 0.000001;
	//int nmax = 1000;

	//res = expansion(ff2T, da, delta_da, alpha, nmax);
	//solution wynik;
	////wynik = fib(ff2T, res[0], res[1], epsilon); 
	//wynik = lag(ff2T, res[0], res[1], epsilon, gamma, nmax);
	//cout << wynik;
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
