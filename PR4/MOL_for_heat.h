#include <bits/stdc++.h>
#include "Function.h"
#include "Matrix.h"
#include "MOL.h"
using namespace std;

class Theta_Method : public MOL_for_heat {
protected:
	double theta;
	virtual void step(const Function<double>& f, const int& n) {
		Matrix<double> coef(m+1, m+1);
		Colvec<double> rhs(m+1);
		for (int i = 1; i < m; ++ i) {
			coef[i][i-1] = -r * theta;
			coef[i][i] = 1 + 2 * theta * r;
			coef[i][i+1] = -r * theta;
			rhs[i] = 	(1-theta) * r * u[n][i-1]
						+ (1 - 2 * (1-theta) * r) * u[n][i]
						+ (1-theta) * r * u[n][i+1];
		}
		coef[0][0] = 1, rhs[0] = 0;
		coef[m][m] = 1, rhs[m] = 0;
		u[n+1] = Gauss_Improved_Solve(coef, rhs);
	}
public:
	Theta_Method(double theta = 0.5) : theta(theta) {}
};

class MOL_heat_with_RK : public MOL_for_heat {
protected:
	int s;
	vector<vector<double>> a;
	vector<double> b;
	int id(int i, int j) {
		return j*(m+1)+i;
	}
	virtual void step(const Function<double>& f, const int& n) {
		Matrix<double> coef(s*(m+1), s*(m+1));
		Colvec<double> rhs(s*(m+1));
		for (int i = 1; i < m; ++ i) {
			for (int j = 0; j < s; ++ j) {
				for (int l = 0; l < s; ++ l) {
					coef[id(i,j)][id(i,l)] += 2*r*a[j][l];
					coef[id(i,j)][id(i-1,l)] -= r*a[j][l];
					coef[id(i,j)][id(i+1,l)] -= r*a[j][l];
				}
				coef[id(i,j)][id(i,j)] += 1;
				rhs[id(i,j)] = r * (u[n][i-1] - 2*u[n][i] + u[n][i+1]) / k;
			}
		}
		for (int j = 0; j < s; ++ j) {
			coef[id(0,j)][id(0,j)] = 1, rhs[id(0,j)] = 0;
			coef[id(m,j)][id(m,j)] = 1, rhs[id(m,j)] = 0;
		}
		Colvec<double> ys = Gauss_Improved_Solve(coef, rhs);
		vector<Colvec<double>> y(s);
		for (int j = 0; j < s; ++ j) y[j] = split(ys, j*(m+1), (j+1)*(m+1));
		u[n+1] = u[n];
		for (int j = 0; j < s; ++ j) u[n+1] += k * b[j] * y[j];
	}
};

class Collocation_Method : public MOL_heat_with_RK {
public:
	Collocation_Method() {
		s = 2;
		a = {
			{5./12, -1./12},
			{3./4, 1./4}
		};
		b = {3./4, 1./4};
	}
};

class Gauss_Method : public MOL_heat_with_RK {
public:
	Gauss_Method() {
		s = 1;
		a = {
			{1./2}
		};
		b = {1};
	}
};