#include <bits/stdc++.h>
#include "function.h"
#include "IVPSolver.h"
using namespace std;

class RK : public IVPSolver{
protected:
	int s, p;
	vector<vector<double>> a;
	vector<double> b, c;
	virtual Colvec<double> step(const Function_nd<double>& f, const Colvec<double>& u, double k) = 0;
public:
	virtual vector<Colvec<double>> Solve(const Function_nd<double>& f, const Colvec<double>& u0, double T, int n, int s = 0, bool omit = 0) {
		long long start = chrono::steady_clock::now().time_since_epoch().count();
		double k = T/n;
		u.resize(n+1);
		u[0] = u0;
		for (int i = 0; i < n; ++ i) u[i+1] = step(f, u[i], k);
		long long end = chrono::steady_clock::now().time_since_epoch().count();
		if (!omit) {
			cout << "============RK Method Return============" << endl;
			cout << "Method: " << method_name << endl;
			cout << "Stages: " << s << endl;
			cout << "Start time: 0, Terminal time: " << T << endl;
			cout << "Steps: " << n << endl;
			cout << "CPU Time: " << (end - start) * 1e-6 << endl;
			cout << "Error: " << ~(u[n] - u[0]) << endl;
			cout << "2-norm of Error: " << vert_2(u[n] - u[0]) << endl;
			Colvec<double> u_half, err(f.n);
			for (int i = 0; i < n; ++ i) {
				u_half = step(f, step(f, u[i], k/2), k/2);
				for (int j = 0; j < f.n; ++ j) err[j] += fabs(u[i+1][j] - u_half[j]) * (1<<p) / ((1<<p)-1);
			}
			cout << "Error Estimated by Richardson extrapolation: " << ~err << endl;
			cout << "2-norm: " << vert_2(err) << endl;
			cout << "=================================================" << endl;
		}
		return u;
	}
};

class ERK : public RK {
protected:
	virtual Colvec<double> step(const Function_nd<double>& f, const Colvec<double>& u, double k) {
		vector<Colvec<double>> y(s);
		for (int i = 0; i < s; ++ i) {
			Colvec<double> r = u;
			for (int j = 0; j < i; ++ j)
				r += k * a[i][j] * y[j];
			y[i] = f(r);
		}
		Colvec<double> res = u;
		for (int i = 0; i < s; ++ i) 
			res += k * b[i] * y[i];
		return res;
	}
};

class IRK : public RK {
protected:
	virtual Colvec<double> step(const Function_nd<double>& f, const Colvec<double>& u, double k) {
		vector<Colvec<double>> y(s);
		for (int i = 0; i < s; ++ i) y[i] = f(u);
		Colvec<double> res = u;
		Colvec<double> last, now;
		do {
			vector<Colvec<double>> ny(s);
			for (int i = 0; i < s; ++ i) {
				Colvec<double> r = u;
				for (int j = 0; j < s; ++ j)
					r += k * a[i][j] * y[j];
				ny[i] = f(r);
			}
			now = u;
			for (int i = 0; i < s; ++ i)
				now += k * b[i] * ny[i];
			last = res;
			res = now;
			y = ny;
		} while(vert_2(last - now) > 1e-14);
		return res;
	}
};

class Classical_4th_RK : public ERK {
public:
	Classical_4th_RK() {
		method_name = "Classical RK";
		s = 4, p = 4;
		a = {
			{0, 0, 0, 0},
			{0.5, 0, 0, 0},
			{0, 0.5, 0, 0},
			{0, 0, 1, 0}
		};
		b = {1/6.0, 1/3.0, 1/3.0, 1/6.0};
		c = {0, 0.5, 0.5, 1};
	}
};

class ESDIRK : public IRK {
public:
	ESDIRK() {
		method_name = "ESDIRK";
		s = 6, p = 4;
		a = {
			{0, 0, 0, 0, 0, 0},
			{1/4.0, 1/4.0, 0, 0, 0, 0},
			{8611/62500.0, -1743/31250.0, 1/4.0, 0, 0, 0},
			{5012029.0/34652500.0, -654441.0/2922500.0, 174375.0/388108.0, 1/4.0, 0, 0},
			{15267082809.0/155376265600.0, -71443401.0/120774400.0, 730878875.0/902184768.0, 2285395.0/8070912.0, 1/4.0, 0},
			{82889/524892.0, 0, 15625/83664.0, 69875/102672.0, -2260/8211.0, 1/4.0}
		};
		b = a[5];
		c = {0, 1/2.0, 83/250.0, 31/50.0, 17/20.0, 1};
	}
};

class Gauss : public IRK {
public:
	Gauss(int s) {
		method_name = "Gauss-Legendre";
		RK::s = s, p = 2*s;
		if (s == 1) {
			a = {{0.5}};
			b = {1};
			c = {0.5};
		}
		else if (s == 2) {
			double q = sqrt(3);
			a = {
				{1.0/4, (3-2*q)/12},
				{(3+2*q)/12, 1.0/4}
			};
			b = {0.5, 0.5};
			c = {(3-q)/6, (3+q)/6};
		}
		else if (s == 3) {
			double q = sqrt(15);
			a = {
				{5.0/36, 2.0/9-q/15, 5.0/36-q/30},
				{5.0/36+q/24, 2.0/9, 5.0/36-q/24},
				{5.0/36+q/30, 2.0/9+q/15, 5.0/36}
			};
			b = {5.0/18, 4.0/9, 5.0/18};
			c = {(5-q)/10, 1.0/2, (5+q)/10};
		}
	}
};