#include<bits/stdc++.h>
#include "function.h"
#include "Dirichlet_Normal.h"
#include "Neuman_Normal.h"
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
	double operator ()(const double& x, const double& y) const {
		return -(u.partial(x, y, 2, 0) + u.partial(x, y, 0, 2));
	}
} f;

class G : public Function_2D<double> {
public:
	double operator ()(const double& x, const double& y) const {
		if (x == 0) return u.partial(0, y, 1, 0);
		if (x == 1) return -u.partial(1, y, 1, 0);
		if (y == 0) return u.partial(x, 0, 0, 1);
		if (y == 1) return -u.partial(x, 1, 0, 1);
		throw 0;
	}
} g;


int main(int argc, char** argv){
	string cmd = argv[1];
	int n = atoi(argv[2]);
	Dirichlet_Normal_Solver D(n, f, u);
	Neuman_Normal_Solver N(n, f, g);
	BVPSolver* A;
	if (cmd == "Dirichlet") A = &D;
	else if (cmd == "Neuman") A = &N;
	A->Solve();
	for (int i = 0; i <= n+1; ++ i, cout << '\n')
		for (int j = 0; j <= n+1; ++ j)
			cout << A->getval(i,j) << ' ';
	double h = 1.0/(n+1);
	for (int i = 0; i <= n+1; ++ i, cout << '\n')
		for (int j = 0; j <= n+1; ++ j)
			cout << u(i*h, j*h) << ' ';
}