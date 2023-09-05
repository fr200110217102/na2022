#include <bits/stdc++.h>
#include "Matrix.h"
#include "Function.h"
#include "MOLFactory.h"
#include "MOL_for_heat.h"
using namespace std;

const int K = 300;
const double pi = acos(-1);

class F : public Function<double>{
public:
	double operator()(const double& x) const {
		if (0.45 < x && x <= 0.5) return 20 * x - 9;
		if (0.5 <= x && x < 0.55) return 11 - 20 * x;
		return 0;
	}
};

double A[K+1];

class U : public Function_2D<double>{
public:
	double operator()(const double& t, const double& x) const {
		double res = 0;
		for (int k = 1; k <= K; ++ k) res += A[k] * exp(-k*k*pi*pi*t) * sin(k*pi*x);
		return res;
	}
} u_rel;

int main(int argc, char** argv) {
	for (int k = 1; k <= K; ++ k) A[k] = 40 * (-sin(0.45*k*pi) + 2*sin(0.5*k*pi) - sin(0.55*k*pi)) / (k*k*pi*pi);

	MOLFactory& Factory = MOLFactory::CreateFactory();
	Factory.Register("FTCS", []() -> unique_ptr<MOL>{return make_unique<Theta_Method>(0);});
	Factory.Register("BTCS", []() -> unique_ptr<MOL>{return make_unique<Theta_Method>(1);});
	Factory.Register("Crank-Nicolson", []() -> unique_ptr<MOL>{return make_unique<Theta_Method>(0.5);});
	Factory.Register("Collocation-Method", []() -> unique_ptr<MOL>{return make_unique<Collocation_Method>();});
	Factory.Register("Gauss-Method", []() -> unique_ptr<MOL>{return make_unique<Gauss_Method>();});

	string name = argv[1];
	double T = atof(argv[2]);
	int n = atoi(argv[3]);
	int m = atoi(argv[4]);

	F f;
	unique_ptr<MOL> Solver = Factory.create(name);

	long long start = chrono::steady_clock::now().time_since_epoch().count();
	vector<Colvec<double>> u = Solver->Solve(f, 1, T, n, m);
	long long end = chrono::steady_clock::now().time_since_epoch().count();

	cout << "====================MOL Returned======================" << endl;
	cout << "Method: " << name << endl;
	cout << "Terminal Time: " << T << endl;
	cout << "Step: " << n << endl;
	cout << "Grid Size: " << m << endl;
	cout << "CPU Time: " << (end - start) * 1e-6 << endl;
	cout << "Terminal Error in 2-norm: ";
	double res = 0;
	for (int i = 1; i < m; ++ i) {
		double e = u[n][i] - u_rel(T, 1.0*i/m);
		res += e*e;
	}
	res = sqrt(res / (m-1));
	cout << res << endl;
	cout << "======================================================" << endl;

	if (argc > 5) {
		ofstream out(argv[5], ios::out);
		out << "x,u\n";
		for (int i = 0; i <= m; ++ i) out << 1.0*i/m << ',' << u[n][i] << '\n';
		out.flush();
		out.close();
	}
}