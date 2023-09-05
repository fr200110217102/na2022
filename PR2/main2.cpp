#include<bits/stdc++.h>
#include "Function.h"
#include "multigrid.h"
using namespace std;

class U : public Function_2D<double> {
public:
	double operator ()(const double& x, const double& y) const {
		return exp(y + sin(x));
	}
	double partial(const double& x, const double& y, const int& i, const int& j) const {
		if (i == 0) return exp(y + sin(x));
		if (i == 1 && j == 0) return cos(x) * exp(y + sin(x));
		if (i == 2 && j == 0) return (cos(x) * cos(x) - sin(x)) * exp(y + sin(x));
		throw 0;
	}
} u;

class F : public Function_2D<double> {
public:
	Function_2D<double>& u;
	F(Function_2D<double>& u) : u(u) {}
	double operator ()(const double& x, const double& y) const {
		return -(u.partial(x, y, 2, 0) + u.partial(x, y, 0, 2));
	}
};

typedef unsigned int Cond_t;
template<Cond_t Cond_type>
class G : public Function_2D<double> {
public:
	Function_2D<double>& u;
	G(Function_2D<double>& u) : u(u) {}
	double operator ()(const double& x, const double& y) const {
		if (y == 0) return Cond_type & 1 ? u.partial(x, 0, 0, 1) : u(x, 0);
		if (x == 0) return Cond_type & 2 ? u.partial(0, y, 1, 0) : u(0, y) ;
		if (y == 1) return Cond_type & 4 ? -u.partial(x, 1, 0, 1) : u(x, 1);
		if (x == 1) return Cond_type & 8 ? -u.partial(1, y, 1, 0) : u(1, y);
		throw 0;
	}
};

int main(int argc, char** argv) {
	int n = atoi(argv[1]);
	const Restriction_method i1 = string(argv[2]) == "I" ? Injection : Full_Weighting;
	const Interpolation_method i2 = string(argv[3]) == "L" ? Linear : Quadratic;
	const Cycle_method i3 = string(argv[4]) == "V" ? V_cycle : FMG;
	int T1 = atoi(argv[5]), T2 = atoi(argv[6]);
	auto Solver_D = Multigrid<2, 0>(F(u), G<0>(u));
	Solver_D.Solve(n, i1, i2, i3, T1, T2, 1e-8, 1, u);
	auto Solver_N = Multigrid<2, 15>(F(u), G<15>(u));
	Solver_N.Solve(n, i1, i2, i3, T1, T2, 1e-8, 1, u);
	auto Solver_M = Multigrid<2, 3>(F(u), G<3>(u));
	Solver_M.Solve(n, i1, i2, i3, T1, T2, 1e-8, 1, u);
}