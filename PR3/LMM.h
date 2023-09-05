#include <bits/stdc++.h>
#include "Function.h"
#include "IVPSolver.h"
using namespace std;

class LMM : public IVPSolver{
protected:
	virtual void step(const Function_nd<double>& f, int n, double k, int s) = 0;
public:
	virtual vector<Colvec<double>> Solve(const Function_nd<double>& f, const Colvec<double>& u0, double T, int n, int s = 0, bool omit = 0) {
		long long start = chrono::steady_clock::now().time_since_epoch().count();
		double k = T/n;
		u.resize(1);
		u[0] = u0;
		if (s>1) {
			if (s>2) {
				double _k = pow(k, 1.0*s/(s-1));
				int m = k / _k;
				vector<Colvec<double>> v = Solve(f, u0, k*(s-1), m*(s-1), s-1, 1);
				u.resize(s);
				for (int i = 1; i <= s-1; ++ i) u[i] = v[i*m];
			}
			else {
				double m = 1;
				while (m < 1/k) m *= 2;
				double kk = k/m;
				u.resize(3);
				step(f, 0, kk, 1);
				for (double i = 1; i < m; i *= 2) {
					step(f, 0, kk*i, 2);
					u[1] = u[2];
				}
			}
		}
		u.resize(n+1);
		for (int i = s; i <= n; ++ i) step(f, i-s, k, s);
		long long end = chrono::steady_clock::now().time_since_epoch().count();
		if (!omit) {
			cout << "============LMM Method Return============" << endl;
			cout << "Method: " << method_name << endl;
			cout << "Stage: " << s << endl;
			cout << "Start time: 0, Terminal time: " << T << endl;
			cout << "Steps: " << n << endl;
			cout << "CPU Time: " << (end - start) * 1e-6 << endl;
			cout << "Error: " << ~(u[n] - u[0]) << endl;
			cout << "2-norm of Error: " << vert_2(u[n] - u[0]) << endl;
			cout << "=========================================" << endl;
		}
		return u;
	}
};

class Adam_Bashforth : public LMM {
protected:
	vector<double> beta[5];
	virtual void step(const Function_nd<double>& f, int n, double k, int s) {
		u[n+s] = u[n+s-1];
		for (int i = 0; i <= s-1; ++ i)
			u[n+s] += k * beta[s][i] * f(u[n+i]);
	}
public:
	Adam_Bashforth()  {
		method_name = "Adam-Bashforth";
		beta[1] = {1};
		beta[2] = {-1/2.0, 3/2.0};
		beta[3] = {5/12.0, -16/12.0, 23/12.0};
		beta[4] = {-9/24.0, 37/24.0, -59/24.0, 55/24.0};
	}
};


class Adam_Moulton : public LMM {
protected:
	vector<double> beta[5];
	virtual void step(const Function_nd<double>& f, int n, double k, int s) {
		u[n+s] = u[n+s-1];
		Colvec<double> last, now;
		do {
			now = u[n+s-1];
			for (int i = 0; i <= s; ++ i)
				now += k * beta[s][i] * f(u[n+i]);
			last = u[n+s];
			u[n+s] = now;
		} while (vert_2(last - now) > 1e-14);
	}
public:
	Adam_Moulton() {
		method_name = "Adam-Moulton";
		beta[1] = {1/2.0, 1/2.0};
		beta[2] = {-1/12.0, 8/12.0, -5/12.0};
		beta[3] = {1/24.0, -5/24.0, 19/24.0, 9/24.0};
		beta[4] = {-19/720.0, 106/720.0, -264/720.0, 646/720.0, 251/720.0};
	}
};

class BDF : public LMM {
protected:
	vector<double> alpha[5];
	double beta[5];
	virtual void step(const Function_nd<double>& f, int n, double k, int s) {
		u[n+s] = u[n+s-1];
		Colvec<double> last, now;
		do {
			now = k * beta[s] * f(u[n+s]);
			for (int i = 0; i <= s-1; ++ i)
				now -= alpha[s][i] * u[n+i];
			last = u[n+s];
			u[n+s] = now;
		} while (vert_2(last - now) > 1e-14);
	}
public:
	BDF() {
		method_name = "BDF";
		beta[1] = 1,		alpha[1] = {-1};
		beta[2] = 2/3.0,	alpha[2] = {1/3.0, -4/3.0};
		beta[3] = 6/11.0,	alpha[3] = {-2/11.0, 9/11.0, -18/11.0};
		beta[4] = 12/25.0,	alpha[4] = {3/25.0, -16/25.0, 36/25.0, -48/25.0};
	}
};