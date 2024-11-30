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
#include"user_funs.h"

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
		#define TEORETYCZNE4

		lab0();
		lab1();
		lab2();
		lab3();
		lab4();
		lab5();
		lab6();
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
	double epsilon = 1e-2;
	int Nmax = 10000;
	matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);

#ifdef TEORETYCZNE0
	
	//Funkcja testowa
	solution opt;
	a(0) = -1;
	a(1) = 2;

	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
	cout << opt << endl << endl;
	solution::clear_calls();

#endif // TEORETYCZNE0

#ifdef PRAKTYCZNE0

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

#endif // PRAKTYCZNE0
}

void lab1()
{
#ifdef TEORETYCZNE1

	std::ofstream Sout("symulacja_lab1.csv");
	
	//zadanie teoretyczne

	double* res = new double[2] { 0, 0 };
	double x0 = 50, d = 5, alpha = 1.5;
	int Nmax = 10000;

	//double a = 50, b = 70;
	double epsilon = 0.0001;
	double gamma = 0.000001;
	solution wynik;

	for (int i = 0; i < 1; i++)
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

#endif // TEORETYCZNE1

#ifdef PRAKTYCZNE1

	//zadanie praktyczne

	double* res = new double[2] { 0, 0 };
	double da = 0.005, delta_da = 0.002;
	double alpha = 1.5, epsilon = 0.0001, gamma = 0.000001;
	int nmax = 1000;

	res = expansion(ff1R, da, delta_da, alpha, nmax);
	solution wynik;
	//wynik = fib(ff1R, res[0], res[1], epsilon); 
	cout << "Metoda Fibonacciego: " << endl;
	cout << "Optymalna wielosc otwou D_A: " << wynik.x << "\nMaksymalna temperatura wody w zbiorniku: " << wynik.y + 50 << "|" << wynik.y << "\nLiczna wywolan fukcji: " << wynik.f_calls << "\nExit flag: " << wynik.flag << endl;;
	//cout << wynik;
	wynik = lag(ff1R, res[0], res[1], epsilon, gamma, nmax);
	cout << "Metoda Lagrangea: " << endl;
	cout << "Optymalna wielosc otworu D_A: " << wynik.x << "\nMaksymalna temperatura wody w zbiorniku: " << wynik.y + 50 << "|" << wynik.y << "\nLiczna wywolan fukcji: " << wynik.f_calls << "\nExit flag: " << wynik.flag << endl;
	//cout << wynik;

#endif // PRAKTYCZNE1

#ifdef SYMULACJA1

	//symlacja 

	// Warunki początkowe
	matrix Y0(3, 1);
	Y0(0) = 5.0;   // Początkowa objętość w zbiorniku A (VA)
	Y0(1) = 1.0;   // Początkowa objętość w zbiorniku B (VB)
	Y0(2) = 20.0;  // Początkowa temperatura w zbiorniku B (TB)

	// Czas symulacji
	double t0 = 0.0;            // Początkowy czas
	double tend = 2000.0;       // Końcowy czas symulacji
	double dt = 1.0;            // Krok czasowy

	matrix ud1(1, 1), ud2;
	ud1(0) = m2d(wynik.x);
	
	matrix* S = solve_ode(df1, t0, dt, tend, Y0, ud1, ud2);

	int N = static_cast<int>(floor((tend - t0) / dt) + 1);
	std::vector<double> t_values(N);
	std::vector<double> VA_values(N);
	std::vector<double> VB_values(N);
	std::vector<double> TB_values(N);

	for (int i = 0; i < N; ++i) {
		t_values[i] = S[0](i);
		VA_values[i] = S[1](i, 0); // Pierwsza kolumna to VA
		VB_values[i] = S[1](i, 1); // Druga kolumna to VB
		TB_values[i] = S[1](i, 2); // Trzecia kolumna to TB
	}

	std::ofstream file("symulacja_lab1.csv");
	file << "Czas;VA;VB;TB\n";
	for (int i = 0; i < N; ++i) {
		file << t_values[i] << "x;" << VA_values[i] << "x;" << VB_values[i] << "x;" << TB_values[i] << "x\n";
	}
	file.close();

#endif // SYMULACJA1

}

void lab2()
{

	srand(time(NULL));
	std::ofstream Sout("symulacja_lab2.csv");

	matrix X;
	double step = 0.01, alpha = 0.8, beta = 0.1, epsilon = 0.0001;
	double a, b;
	int Nmax = 1000;

#ifdef TEORETYCZNE2

	// zadanie teoretyczne

	for (int i = 0; i < 100; i++)
	{
		a = ((rand() % 200) / 100.0) - 1;
		b = ((rand() % 200) / 100.0) - 1;
		alpha = 0.8;
		X = matrix(2, new double[2] {a, b});
		solution hooke = HJ(ff2T, X, step, alpha, epsilon, Nmax);
		//Sout << "x" << a << ";" << "x" << b << ";" << "x" << hooke.x(0) << ";" << "x" << hooke.x(1) << ";" << "x" << hooke.y << ";" << "x" << solution::f_calls << ";";
		//cout << hooke;

		
		alpha = 1.8;
		matrix Step = matrix(2, new double[2] { step, step});
		solution rosen = Rosen(ff2T, X, Step, alpha, beta, epsilon, Nmax);
		//Sout << "x" << rosen.x(0) << ";" << "x" << rosen.x(1) << ";" << "x" << rosen.y << ";" << "x" << solution::f_calls << "\n";
		//cout << rosen;
	}
	
#endif //TEORETYCZNE2

#ifdef PRAKTYCZNE2

	//problem rzeczywisty

	//sprawdzenie poprawnosci
	matrix x(2, 1, 5); // macierz 2x1 wypelniona wartoscia 5
	cout << ff2R(x);
	
	X = matrix(2, new double[2] {5, 5});
	solution wynikHJ = HJ(ff2R, X, step, alpha, epsilon, Nmax);
	//cout << wynikHJ;

	alpha = 1.8;
	matrix Step = matrix(2, new double[2] { step, step});
	solution wynikR = Rosen(ff2R, X, Step, alpha, beta, epsilon, Nmax);
	//cout << wynikR;

	//symulacja 

	double t0 = 0.0;
	double tend = 100.0;
	double dt = 0.1;

	// Initial conditions for the state vector Y
	matrix Y0(2, 1); 
	Y0(0) = 0.0;     
	Y0(1) = 0.0;     

	// Reference values for the desired angle and angular velocity
	matrix ud1(2, 1);
	ud1(0) = 3.14;   
	ud1(1) = 0.0;    

	// Gain parameters, k1 and k2, within the range [0, 10]
	matrix ud2(2, 1);
	//ud2(0) = wynikHJ.x(0);     
	ud2(0) = wynikR.x(0);     
	//ud2(1) = wynikHJ.x(1);     
	ud2(1) = wynikR.x(1);     

	matrix* result = solve_ode(df2, t0, dt, tend, Y0, ud1, ud2);

	int n = get_len(result[0]);
	std::cout << "Time\tAngle\tAngular Velocity" << std::endl;
	for (int i = 0; i < n; i++) {
		std::cout << result[0](i) << "\t"      
			<< result[1](i, 0) << "\t"  
			<< result[1](i, 1) << std::endl; 
	}

	std::ofstream file("symulacja_lab2.csv");
	file << "Czas;Kat;Predkosc katowa\n";
	for (int i = 0; i < n; ++i) {
		file << result[0](i) << "x;" << result[1](i, 0) << "x;" << result[1](i, 1) << "x\n";
	}
	file.close();

#endif // PRAKTYCZNE2
}

void lab3()
{
	// Dane dokładnościowe
	double epsilon = 1E-3;
	int Nmax = 10000;
	double c_inside = 100;
	double dc_inside = 0.2;
	double c_outside = 1.0;
	double dc_outside = 1.5;

	std::ofstream Sout("symulacja_lab3.csv");

#ifdef TEORETYCZNE3

	// Nagłówki w pliku CSV
	Sout << "x0_1;x0_2;x1_out;x2_out;norm_out;y_out;f_calls_out;x1_in;x2_in;norm_in;y_in;f_calls_in\n";

	// Generator liczb losowych
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> x0_dist(1.5, 5.5);

	// Stringstream do zapisu danych
	std::stringstream test_ss;

	// Rozwiązanie dla wyników testowych
	solution test_sol;

	// Dane a dla testów
	matrix a = matrix(4.0);

	// Punkty startowe dla testów
	matrix test_x0{};

	for (int i = 0; i < 3; ++i) {
		if (i == 0)
			a = matrix(4.0);
		else if (i == 1)
			a = matrix(4.4934);
		else
			a = matrix(5.0);

		for (int j = 0; j < 100; ++j) {
			test_x0 = matrix(2, new double[2] {x0_dist(gen), x0_dist(gen)});
			test_ss << test_x0(0) << "X ;" << test_x0(1) << "X ;";

			// Zewnętrzne rozwiązanie
			test_sol = pen(ff3T_outside, test_x0, c_outside, dc_outside, epsilon, Nmax, a);
			//cout << test_sol;
			test_ss << test_sol.x(0) << "X ;" << test_sol.x(1) << "X ;"
				<< sqrt(pow(test_sol.x(0), 2) + pow(test_sol.x(1), 2)) << "X ;"
				<< test_sol.y << "X ;" << test_sol.f_calls << "X ;";
			solution::clear_calls();

			// Wewnętrzne rozwiązanie
			test_sol = pen(ff3T_inside, test_x0, c_inside, dc_inside, epsilon, Nmax, a);
			//cout << test_sol;
			test_ss << test_sol.x(0) << "X ;" << test_sol.x(1) << "X ;"
				<< sqrt(pow(test_sol.x(0), 2) + pow(test_sol.x(1), 2)) << "X ;"
				<< test_sol.y << "X ;" << test_sol.f_calls << "X \n";
			solution::clear_calls();
		}
	}

	// Zapis zawartości stringstream do pliku
	Sout << test_ss.str();
	Sout.close();

#endif // TEORETYCZNE3

#ifdef PRAKTYCZNE3

	// Nagłówki w pliku CSV
	Sout << "time;x_position;y_position\n";

	//Dane zadania
	matrix ud1 = matrix(5, new double[5] {
		0.47, //Współczynnik oporu (C) [-]
		1.2, //Gęstość powietrza (rho) [kg/m^3]
		0.12, //Promień piłki (r) [m]
		0.6, //Masa piłki (m) [kg]
		9.81 //Przyśpieszenie ziemskie (g) [m/s^2]
});

	//Początkowe wartości szukania minimum
	matrix x0 = matrix(2, new double[2] {-4.0, 2.0});

	//Szukanie optymalnej prędkości początkowej po osi x i początkowej prędkości obrotowej
	solution opt = pen(ff3R, x0, c_outside, dc_outside, epsilon, Nmax, ud1);
	std::cout << opt << "\n";

	//Symulacja lotu piłki dla wyznaczonych ograniczeń
	matrix Y0(4, new double[4] {0.0, opt.x(0), 100, 0});
	matrix* Y = solve_ode(df3, 0.0, 0.01, 7.0, Y0, ud1, opt.x(1));

	// Zapis wyników symulacji do pliku CSV
	int num_rows = get_len(Y[0]);
	for (int i = 0; i < num_rows; ++i) {
		Sout << Y[0](i) << "x ;"  // Czas
			<< Y[1](i, 0) << "x ;"  // Pozycja w osi x
			<< Y[1](i, 2) << "x \n";  // Pozycja w osi y
	}

	// Zamknięcie pliku i usunięcie dynamicznie alokowanej pamięci
	Sout.close();
	delete[] Y;
	

#endif // PRAKTYCZNE3

}

void lab4()
{
#ifdef TEORETYCZNE4
	double epsilon = 1e-6;
	double h0 = 0.05;
	double a = 0.0;
	double b = 1.0;
	int Nmax = 1000;

	matrix x0 = matrix(2, new double[2] {0.0, 0.0});
	matrix ud1 = NAN;
	matrix ud2 = NAN;

	//solution grad_result = SD(ff4T, gf4T, x0, h0, epsilon, Nmax, ud1, ud2);
	//solution grad_result = CG(ff4T, gf4T, x0, h0, epsilon, Nmax, ud1, ud2);
	solution grad_result = Newton(ff4T, gf4T, hf4T, x0, h0, epsilon, Nmax, ud1, ud2);
	std::cout << grad_result << "\n";
#endif
}

void lab5()
{

}

void lab6()
{

}
