#include <bits/stdc++.h>
#include "Matrix.h"
#include "Function.h"
#include "MOLFactory.h"
#include "MOL_for_adv.h"
using namespace std;

const int K = 300;
const double pi = acos(-1);

class F : public Function<double>{
public:
	double operator()(const double& x) const {
		return exp(-20*(x-2)*(x-2)) + exp(-(x-5)*(x-5));
	}
};

class U : public Function_2D<double>{
public:
	double operator()(const double& t, const double& x) const {
		return exp(-20*(x-t-2)*(x-t-2)) + exp(-(x-t-5)*(x-t-5));
	}
} u_rel;

int main(int argc, char** argv) {
	MOLFactory& Factory = MOLFactory::CreateFactory();
	Factory.Register("LeapFrog", []() -> unique_ptr<MOL>{return make_unique<LeapFrog>();});
	Factory.Register("Lax-Friedrichs", []() -> unique_ptr<MOL>{return make_unique<Lax_Friedrichs>();});
	Factory.Register("Lax-Wendroff", []() -> unique_ptr<MOL>{return make_unique<Lax_Wendroff>();});
	Factory.Register("Upwind", []() -> unique_ptr<MOL>{return make_unique<Upwind>();});
	Factory.Register("Beam-Warming", []() -> unique_ptr<MOL>{return make_unique<Beam_Warming>();});

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
		double e = u[n][i] - u_rel(T, 25.*i/m);
		res += e*e;
	}
	res = sqrt(25. * res / (m-1));
	cout << res << endl;
	cout << "======================================================" << endl;

	if (argc > 5) {
		ofstream out(argv[5], ios::out);
		out << "x,u\n";
		for (int i = 0; i <= m; ++ i) out << 25.0*i/m << ',' << u[n][i] << '\n';
		out.flush();
		out.close();
	}
}