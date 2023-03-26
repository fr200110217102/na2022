#include<bits/stdc++.h>
#include "function.h"
#include "geo_2D.h"
#include "Normal_BVPSolver.h"
#include "Irnormal_BVPSolver.h"
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

class V : public Function_2D<double> {
public:
	double operator ()(const double& x, const double& y) const {
		return exp(-(x*x+y*y));
	}
	double partial(const double& x, const double& y, const int& i, const int& j) const {
		if (i == 1 && j == 0) return -2 * x * exp(-(x*x+y*y));
		if (i == 2 && j == 0) return (4*x*x - 2) * exp(-(x*x+y*y));
		if (i == 0 && j == 1) return -2 * x * exp(-(x*x+y*y));
		if (i == 0 && j == 2) return (4*y*y - 2) * exp(-(x*x+y*y));
		throw 0;
	}
} v;

class W : public Function_2D<double> {
public:
	double operator ()(const double& x, const double& y) const {
		return 1/sqrt(x*x+y*y);
	}
	double partial(const double& x, const double& y, const int& i, const int& j) const {
		if (i == 1 && j == 0) return -x / (x*x+y*y) / sqrt(x*x+y*y);
		if (i == 2 && j == 0) return (2*x*x-y*y) / (x*x+y*y) / (x*x+y*y) / sqrt(x*x+y*y);
		if (i == 0 && j == 1) return -y / (x*x+y*y) / sqrt(x*x+y*y);
		if (i == 0 && j == 2) return (2*y*y-x*x) / (x*x+y*y) / (x*x+y*y) / sqrt(x*x+y*y);
		throw 0;
	}
} w;

class F : public Function_2D<double> {
public:
	Function_2D<double>& u;
	F(Function_2D<double>& u) : u(u) {}
	double operator ()(const double& x, const double& y) const {
		return -(u.partial(x, y, 2, 0) + u.partial(x, y, 0, 2));
	}
};

class G : public Function_2D<double> {
public:
	Function_2D<double>& u;
	string cond;
	G(Function_2D<double>& u, const string& s) : u(u), cond(s) {}
	double operator ()(const double& x, const double& y) const {
		if (y == 0) return cond[0] == 'D' ? u(x, 0) : u.partial(x, 0, 0, 1);
		if (x == 0) return cond[1] == 'D' ? u(0, y) : u.partial(0, y, 1, 0);
		if (y == 1) return cond[2] == 'D' ? u(x, 1) : -u.partial(x, 1, 0, 1);
		if (x == 1) return cond[3] == 'D' ? u(1, y) : -u.partial(1, y, 1, 0);
		throw 0;
	}
};

class G1 : public Function_2D<double> {
public:
	Function_2D<double>& u;
	circle D;
	string cond;
	G1(Function_2D<double>& u, double x0, double y0, double r, const string& s) : u(u), D(x0, y0, r), cond(s) {}
	double operator ()(const double& x, const double& y) const {
		if (y == 0) return cond[0] == 'D' ? u(x, 0) : u.partial(x, 0, 0, 1);
		if (x == 0) return cond[1] == 'D' ? u(0, y) : u.partial(0, y, 1, 0);
		if (y == 1) return cond[2] == 'D' ? u(x, 1) : -u.partial(x, 1, 0, 1);
		if (x == 1) return cond[3] == 'D' ? u(1, y) : -u.partial(1, y, 1, 0);
		if (pnt_to_circle(pnt(x,y), D) == 0) {
			if (cond[4] == 'D') return u(x, y);
			double arg = D.arg(pnt(x,y));
			return cos(arg) * u.partial(x, y, 1, 0) + sin(arg) * u.partial(x, y, 0, 1);
		}
		throw 0;
	}
};

int main(int argc, char** argv){
	string cmd0 = argv[1], cmd1 = argv[2], cmd2 = argv[3];
	Function_2D<double>& t = u;
	if (cmd0 == "v") t = v;
	if (cmd0 == "w") t = w;
	int n = atoi(argv[4]);
	F f(t);
	if (cmd1 == "Normal") {
		if (cmd2 == "Dirichlet") cmd2 = "DDDD";
		if (cmd2 == "Neumann") cmd2 = "NNNN";
		G g(t, cmd2);
		Normal_BVPSolver M(n, f, g, cmd2);
		M.Solve();
		M.Summary(t);
	}
	else {
		double x0 = atof(argv[5]);
		double y0 = atof(argv[6]);
		double r = atof(argv[7]);
		if (cmd2 == "Dirichlet") cmd2 = "DDDDD";
		if (cmd2 == "Neumann") cmd2 = "NNNNN";
		G1 g(t, x0, y0, r, cmd2);
		Irnormal_BVPSolver MI(n, f, g, x0, y0, r, cmd2);
		MI.Solve();
		MI.Summary(t);
	}
}