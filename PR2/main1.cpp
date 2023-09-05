#include<bits/stdc++.h>
#include "Function.h"
#include "multigrid.h"
using namespace std;

class U : public Function<double> {
public:
	virtual double operator ()(const double& x) const {
		return exp(sin(x));
	}
	virtual double d(const double& x, const int& k = 1) const {
		if (k == 1) return cos(x) * exp(sin(x));
		if (k == 2) return (-sin(x) + cos(x)*cos(x)) * exp(sin(x));
		throw 0;
	}
}u;

class F : public Function<double> {
public:
	Function<double>& u;
	F(Function<double>& u) : u(u) {}
	virtual double operator()(const double& x) const {
		return -u.d(x, 2);
	}
};

typedef unsigned int Cond_t;
template<Cond_t Cond_type>
class G : public Function<double> {
public:
	Function<double>& u;
	G(Function<double>& u) : u(u) {}
	virtual double operator()(const double& x) const {
		if (x == 0){
			if (Cond_type & 1) return u.d(0);
			else return u(0);
		}
		if (x == 1){
			if (Cond_type & 2) return -u.d(1);
			else return u(1);
		}
		cerr << "Undefined!" << endl;
		exit(-1);
	}
};

int main(int argc, char** argv){
	int n = atoi(argv[1]);
	const Restriction_method i1 = string(argv[2]) == "I" ? Injection : Full_Weighting;
	const Interpolation_method i2 = string(argv[3]) == "L" ? Linear : Quadratic;
	const Cycle_method i3 = string(argv[4]) == "V" ? V_cycle : FMG;
	int T1 = atoi(argv[5]), T2 = atoi(argv[6]);
	auto Solver_D = Multigrid<1, 0>(F(u), G<0>(u));
	Solver_D.Solve(n, i1, i2, i3, T1, T2, 1e-8, 1, u);
	auto Solver_N = Multigrid<1, 3>(F(u), G<3>(u));
	Solver_N.Solve(n, i1, i2, i3, T1, T2, 1e-8, 1, u);
	auto Solver_M = Multigrid<1, 1>(F(u), G<1>(u));
	Solver_M.Solve(n, i1, i2, i3, T1, T2, 1e-8, 1, u);
}