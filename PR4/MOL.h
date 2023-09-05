#ifndef __MOL__
#define __MOL__
#include <bits/stdc++.h>
#include "Function.h"
#include "Matrix.h"

class MOL {
protected:
	int n, m;
	double k, h;
	vector<Colvec<double>> u;
	virtual void init(const Function<double>& f, const double& nu, const double& T, const int& n, const int& m) = 0;
	virtual void step(const Function<double>& f, const int& n) = 0;
public:
	virtual vector<Colvec<double>> Solve(const Function<double>& f, const double& nu, const double& T, const int& n, const int& m) {
		init(f, nu, T, n, m);
		for (int i = 1; i <= n; ++ i) step(f, i-1);
		return u;
	}
};

class MOL_for_heat : public MOL {
protected:
	double r;
	virtual void init(const Function<double>& f, const double& nu, const double& T, const int& n, const int& m) {
		this->n = n, this->m = m;
		k = T/n, h = 1.0/m, r = nu * k / (h*h);
		u.resize(n+1);
		for (int i = 0; i <= n; ++ i) u[i] = Colvec<double>(m+1);
		for (int i = 0; i <= m; ++ i) u[0][i] = f(i*h);
	}
};

class MOL_for_adv : public MOL {
protected:
	double mu;
	virtual void init(const Function<double>& f, const double& a, const double& T, const int& n, const int& m) {
		this->n = n, this->m = m;
		k = T/n, h = 25./m, mu = a * k / h;
		u.resize(n+1);
		for (int i = 0; i <= n; ++ i) u[i] = Colvec<double>(m+1);
		for (int i = 0; i <= m; ++ i) u[0][i] = f(i*h);
	}
};

#endif