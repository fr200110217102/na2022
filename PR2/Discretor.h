#ifndef DISCRETOR
#define DISCRETOR
#include<bits/stdc++.h>
#include "function.h"
#include "Matrix.h"
#include "Sparsed_Matrix.h"
using namespace std;

/*
	Cond_type:
		Boundary Condition. An Integer.
		The kth bit is 0 if the corresponding boundary condition is Dirichlet and 1 if that is Neumann.
*/

typedef unsigned int Cond_t;
template <int dim, Cond_t Cond_type>
class Discretor {};

template <Cond_t Cond_type>
class Discretor<1, Cond_type> {
private:
	int n;					// The grid number.
	double h;				// The grid size.
	const Function<double>& f;	// The rhs function.
	const Function<double>& g;	// The boundary function.
public:
	Sparsed_Matrix<double> coef;	// The coefficients.
	Colvec<double> rhs;		// The right terms.
private:
	void Normal_Laplace_Discretor_1D(int i0, int i1, int i2) {
		coef[i0][i0] = 2 / (h*h);
		coef[i0][i1] = -1 / (h*h);
		coef[i0][i2] = -1 / (h*h);
		rhs[i0] = f(i0*h);
	}
	void Normal_Dirichlet_Discretor_1D(int i0) {
		coef[i0][i0] = 1;
		rhs[i0] = g(i0*h);
	}
	void Normal_Neumann_Discretor_1D(int i0, int i1, int i2) {
		coef[i0][i0] = -1.5 / h;
		coef[i0][i1] = 2 / h;
		coef[i0][i2] = -0.5 / h;
		rhs[i0] = g(i0*h);
	}
public:
	Discretor() : n(n), f(_0<double>()), g(_0<double>()), coef(0, 0), rhs(0) {}
	Discretor(int n, const Function<double>& f, const Function<double>& g)
		:n(n), h(1.0/n), f(f), g(g), coef(n+1, n+1), rhs(n+1) {
		h = 1.0/n;
		if (Cond_type & 1) Normal_Neumann_Discretor_1D(0, 1, 2);
		else Normal_Dirichlet_Discretor_1D(0);
		if (Cond_type & 2) Normal_Neumann_Discretor_1D(n, n-1, n-2);
		else Normal_Dirichlet_Discretor_1D(n);
		for (int i = 1; i < n; ++ i) Normal_Laplace_Discretor_1D(i, i-1, i+1);
	//	if (Cond_type == 3) coef[n/2].clear(), coef[n/2][n/2] = 1, rhs[n/2] = 0;
	}
};

template <Cond_t Cond_type>
class Discretor<2, Cond_type> {
private:
	int n;					// The grid number.
	double h;				// The grid size.
	const Function_2D<double>& f;	// The rhs function.
	const Function_2D<double>& g;	// The boundary function.
public:
	Sparsed_Matrix<double> coef;	// The coefficients.
	Colvec<double> rhs;		// The right terms.
private:
	int id(int i, int j) {
		return i * (n+1) + j;
	}
	void Normal_Laplace_Discretor_2D(int i0, int i1, int i2, int i3, int i4) {
		coef[i0][i0] += 4 / (h*h);
		coef[i0][i1] += -1 / (h*h);
		coef[i0][i2] += -1 / (h*h);
		coef[i0][i3] += -1 / (h*h);
		coef[i0][i4] += -1 / (h*h);
		rhs[i0] += f(i0/(n+1)*h, i0%(n+1)*h);
	}
	void Normal_Dirichlet_Discretor_2D(int i0) {
		coef[i0][i0] += 1;
		rhs[i0] += g(i0/(n+1)*h, i0%(n+1)*h);
	}
	void Normal_Neumann_Discretor_2D(int i0, int i1, int i2) {
		coef[i0][i0] += -1.5 / h;
		coef[i0][i1] += 2 / h;
		coef[i0][i2] += -0.5 / h;
		rhs[i0] += g(i0/(n+1)*h, i0%(n+1)*h);
	}
public:
	Discretor() : n(0), f(_0_2D<double>()), g(_0_2D<double>()), coef(0, 0), rhs(0) {}
	Discretor(int n, const Function_2D<double>& f, const Function_2D<double>& g)
		:n(n), h(1.0/n), f(f), g(g), coef((n+1)*(n+1), (n+1)*(n+1)), rhs((n+1)*(n+1)) {
		h = 1.0/n;
		for (int i = 0; i <= n; ++ i) {
			if (Cond_type & 1) Normal_Neumann_Discretor_2D(id(i,0), id(i,1), id(i,2));
			else Normal_Dirichlet_Discretor_2D(id(i,0));
		}
		for (int j = 0; j <= n; ++ j) {
			if (Cond_type & 2) Normal_Neumann_Discretor_2D(id(0,j), id(1,j), id(2,j));
			else Normal_Dirichlet_Discretor_2D(id(0,j));
		}
		for (int i = 0; i <= n; ++ i) {
			if (Cond_type & 4) Normal_Neumann_Discretor_2D(id(i,n), id(i,n-1), id(i,n-2));
			else Normal_Dirichlet_Discretor_2D(id(i,n));
		}
		for (int j = 0; j <= n; ++ j) {
			if (Cond_type & 8) Normal_Neumann_Discretor_2D(id(n,j), id(n-1,j), id(n-2,j));
			else Normal_Dirichlet_Discretor_2D(id(n,j));
		}
		for (int i = 1; i < n; ++ i)
			for (int j = 1; j < n; ++ j)
				Normal_Laplace_Discretor_2D(id(i,j), id(i-1,j), id(i+1,j), id(i,j-1), id(i,j+1));
	//	if (Cond_type == 15) coef[id(n/2,n/2)].clear(), coef[id(n/2,n/2)][id(n/2,n/2)] = 1, rhs[id(n/2,n/2)] = 0;
	}
};
#endif