#include <bits/stdc++.h>
#include "Function.h"
#include "IVPSolver.h"
using namespace std;

class Emb_RK : public IVPSolver{
	#define err_tolerance 1e-14
	#define rho_max 5.0
	#define rho_min 0.2
	#define rho 0.8
	#define k_max 1.0
protected:
	int s, p[2];
	vector<vector<double>> a;
	vector<double> b[2], c;
	vector<double> t;
	Colvec<double> step(const Function_nd<double>& f, const Colvec<double>& u, double k, bool o) {
		vector<Colvec<double>> y(s);
		for (int i = 0; i < s; ++ i) {
			Colvec<double> r = u;
			for (int j = 0; j < i; ++ j)
				r += k * a[i][j] * y[j];
			y[i] = f(r);
		}
		Colvec<double> res = u;
		for (int i = 0; i < s; ++ i) 
			res += k * b[o][i] * y[i];
		return res;
	}
	double E_ind(const Colvec<double>& u, const Colvec<double>& e) {
		double res = 0;
		for (int i = 0; i < e.row(); ++ i)
			res += pow(fabs(e[i]) / (err_tolerance * (1 + fabs(u[i]))), 2);
		return sqrt(res / e.row());
	}
public:
	virtual vector<Colvec<double>> Solve(const Function_nd<double>& f, const Colvec<double>& u0, double T, int n = 10000, int s = 0, bool omit = 0) {
		long long start = chrono::steady_clock::now().time_since_epoch().count();
		double k = T/n;
		u.push_back(u0);
		t.push_back(0);
		n = 0;
		while (t[n] < T) {
			if (t[n] + k > T) k = T - t[n];
			Colvec<double> u0, u1;
			double e;
			do{
				u0 = step(f, u[n], k, 0);
				u1 = step(f, u[n], k, 1);
				e = E_ind(u[n], u0-u1);
				k *= min(rho_max, max(rho_min, rho * pow(e, -1.0/(min(p[0], p[1])+1))));
				k = min(k, k_max);
			}while(e > 1);
			u.push_back(u0);
			t.push_back(t[n] + k);
			++n;
		}
		long long end = chrono::steady_clock::now().time_since_epoch().count();
		if (!omit) {
			cout << "============Embedded RK Method Return============" << endl;
			cout << "Method: " << method_name << endl;
			cout << "Stage: " << s << endl;
			cout << "Start time: 0, Terminal time: " << T << endl;
			cout << "Steps: " << n << endl;
			cout << "CPU Time: " << (end - start) * 1e-6 << endl;
			cout << "Error: " << ~(u[n] - u[0]) << endl;
			cout << "2-norm of Error: " << vert_2(u[n] - u[0]) << endl;
			cout << "========================================" << endl;
		}
		return u;
	}
};

class Fehlberg_Emb_RK : public Emb_RK {
public:
	Fehlberg_Emb_RK() {
		method_name = "Fehlberg 4(5) embedded RK";
		s = 6, p[0] = 4, p[1] = 5;
		a = {
			{0, 0, 0, 0, 0, 0},
			{1.0/4, 0, 0, 0, 0, 0},
			{3.0/32, 9.0/32, 0, 0, 0, 0},
			{1932.0/2197, -7200.0/2197, 7296.0/2197, 0, 0, 0},
			{439.0/216, -8, 3680.0/513, -845.0/4104, 0, 0},
			{-8.0/27, 2, -3544.0/2565, 1859.0/4104, -11.0/40, 0},
		};
		b[0] = {25.0/216, 0, 1408.0/2565, 2197.0/4104, -1.0/5, 0};
		b[1] = {16.0/135, 0, 6656.0/12825, 28561.0/56430, -9.0/50, 2.0/55};
		c = {0, 1.0/4, 3.0/8, 12.0/13, 1, 1.0/2};
	}
};

class Dormand_Prince_Emb_RK : public Emb_RK {
public:
	Dormand_Prince_Emb_RK() {
		method_name = "Dormand-Prince 5(4) embedded RK";
		s = 7, p[0] = 5, p[1] = 4;
		a = {
			{0, 0, 0, 0, 0, 0, 0},
			{1.0/5, 0, 0, 0, 0, 0, 0},
			{3.0/40, 9.0/40, 0, 0, 0, 0, 0},
			{44.0/45, -56.0/15, 32.0/9, 0, 0, 0, 0},
			{19372.0/6561, -25360.0/2187, 64448.0/6561, -212.0/729, 0, 0, 0},
			{9017.0/3168, -355.0/33, 46732.0/5247, 49.0/176, -5103.0/18656, 0, 0},
			{35.0/384, 0, 500.0/1113, 125.0/192, -2187.0/6784, 11.0/84, 0}
		};
		b[0] = a[6];
		b[1] = {5179.0/57600, 0, 7571.0/16695, 393.0/640, -92097.0/339200, 187.0/2100, 1.0/40};
		c = {0, 1.0/5, 3.0/10, 4.0/5, 8.0/9, 1, 1};
	}
};