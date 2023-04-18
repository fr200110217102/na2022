#include<bits/stdc++.h>
#include "function.h"
#include "Matrix.h"
#include "Sparsed_Matrix.h"
#include "Discretor.h"
using namespace std;

/*
	Cond:
		Boundary Condition. An Integer.
		The kth bit is 0 if the corresponding boundary condition is Dirichlet and 1 if that is Neumann.
*/
typedef unsigned int Cond_t;

// Restriction: Full Weighting and Injection.
enum Restriction_method {Full_Weighting, Injection};
// Restriction for 1D.
// Input: A vector with n elements.
// Output: A vector with n/2 elements.
Colvec<double> Injection_Restriction_1D(int n, const Colvec<double>& v) {
	Colvec<double> u(n/2+1);
	for (int i = 0; i <= n/2; ++ i)
		u[i] = v[i*2];
	return u;
}
Colvec<double> Full_Weighting_Restriction_1D(int n, const Colvec<double>& v) {
	Colvec<double> u(n/2+1);
	u[0] = v[0], u[n/2] = v[n];
	for (int i = 1; i <= n/2-1; ++ i)
		u[i] = (v[i*2-1] + 2*v[i*2] + v[i*2+1]) / 4;
	return u;
}

// Interpolation: Linear and Quadratic.
enum Interpolation_method {Linear, Quadratic};
// Interpolation for 1D.
// Input: A vector with n elements.
// Output: A vector with n*2 elements.
Colvec<double> Linear_Interpolation_1D(int n, const Colvec<double>& v) {
	Colvec<double> u(n*2+1);
	for (int i = 0; i < n; ++ i) {
		u[i*2] = v[i];
		u[i*2+1] = (v[i] + v[i+1]) / 2;
	}
	u[n*2] = v[n];
	return u;
}
Colvec<double> Quadratic_Interpolation_1D(int n, const Colvec<double>& v) {
	Colvec<double> u(n*2+1);
	for (int i = 0; i < n-1; ++ i) {
		u[i*2] = v[i];
		u[i*2+1] = (3*v[i] + 6*v[i+1] - v[i+2])/ 8;
	}
	u[n*2-2] = v[n-1];
	u[n*2-1] = (-v[n-2] + 6*v[n-1] + 3*v[n]) / 8;
	u[n*2] = v[n];
	return u;
}

enum Cycle_method {V_cycle, FMG};

// Definite the dim and boundary condition type as the template parameters.
template<int dim, Cond_t Cond_type>
class Multigrid {};

template <Cond_t Cond_type>
class Multigrid <1, Cond_type> {
private:
	double w;									// w is the release coefficient.
	const Function<double>& f;					// The rhs function.
	const Function<double>& g;					// The boundary function.
	Restriction_method Restriction_type;
	Interpolation_method Interpolation_type;
	Colvec<double> (*Restriction) (int, const Colvec<double>&);
	Colvec<double> (*Interpolation) (int ,const Colvec<double>&);
	Cycle_method Cycle_type;					// Cycle: V-cycle and FMG.
	map<int, Discretor<1, Cond_type>> D;			// Discretors for different grids.
	Colvec<double> sol;							// The solution vector.
public:
	Multigrid(const Function<double>& f, const Function<double>& g, Restriction_method Restriction_type, Interpolation_method Interpolation_type, Cycle_method Cycle_type):
		w(2.0/3), f(f), g(g), Restriction_type(Restriction_type), Interpolation_type(Interpolation_type), Cycle_type(Cycle_type) {
			if (Restriction_type == Injection) Restriction = Injection_Restriction_1D;
			else if (Restriction_type == Full_Weighting) Restriction = Full_Weighting_Restriction_1D;
			if (Interpolation_type == Linear) Interpolation = Linear_Interpolation_1D;
			else if (Interpolation_type == Quadratic) Interpolation = Quadratic_Interpolation_1D;
		}
private:
	Colvec<double> Jacobi(int n, Sparsed_Matrix<double>& A, Colvec<double>& b, const Colvec<double>& v0, int T) {
		Colvec<double> v(v0);
		Sparsed_Matrix<double> trns(n+1, n+1);
		Colvec<double> c(n+1);
		for (int i = 0; i <= n; ++ i) {
			double dii = A[i][i];
			for (auto & [j, x] : A[i])
				if (j != i) trns[i][j] = -x / dii;
			c[i] = b[i] / dii;
		}
		Colvec<double> u(n+1);
		for (int i = 0; i < T; ++ i) {
			u = trns * v + c;
			v = w * u + (1-w) * v;
		}
		return v;
	}

	// v: Initial guess.
	// A: Discrete matrix for grid number n.
	// b: Discrete rhs for grid number n. 
	Colvec<double> VC(int n, Sparsed_Matrix<double>& A, Colvec<double>& b, Colvec<double>& v0, int T1, int T2) {
		Colvec<double> v = Jacobi(n, A, b, v0, T1);
		if (n <= 2) return Jacobi(n, A, b, v, T2);
		Colvec<double> c = Restriction(n, b - A * v);
		if (!D.count(n/2)) D.insert({n/2, Discretor<1, Cond_type>(n/2, f, g)});
		Colvec<double> zero(n/2+1);
		Colvec<double> v1 = VC(n/2, D[n/2].coef, c, zero, T1, T2);
		v = v + Interpolation(n/2, v1);
		return Jacobi(n, A, b, v, T2);
	}
	Colvec<double> FMGC(int n, Sparsed_Matrix<double>& A, Colvec<double>& b, int T1, int T2) {
		Colvec<double> zero(n+1);
		if (n <= 2) return VC(n, A, b, zero, T1, T2);
		Colvec<double> c = Restriction(n, b);
		if (!D.count(n/2)) D.insert({n/2, Discretor<1, Cond_type>(n/2, f, g)});
		Colvec<double> v = Interpolation(n/2, FMGC(n/2, D[n/2].coef, c, T1, T2));
		return VC(n, A, b, v, T1, T2);
	}
public:
	// n, the finest grid, must be a power of 2.
	// T1 and T2 are the Two-Grid iteration times specified by the user.
	void Solve(int n, int T1 = 5, int T2 = 5, double eps = 1e-8, bool is_test = 0, const Function<double>& real = _0<double>()) {
		if (!D.count(n)) D.insert({n, Discretor<1, Cond_type>(n, f, g)});
		Colvec<double> zero(n+1);
		Colvec<double> rhs = D[n].rhs;
		sol = zero;
		int iter = 0;
		if (Cycle_type == V_cycle) {
			while (1) {
				++ iter;
				Colvec<double> dta = VC(n, D[n].coef, rhs, zero, T1, T2);
				sol = sol + dta;
			//	cout << iter << ' ' << vert_2(D[n].coef * dta) / sqrt(n+1) << endl;
				if (vert_2(D[n].coef * dta) / sqrt(n+1) < eps) break;
				rhs = rhs - D[n].coef * dta;
			}
		}
		else {
			while (1) {
				++ iter;
				Colvec<double> dta = FMGC(n, D[n].coef, rhs, T1, T2);
				sol = sol + dta;
			//	cout << iter << ' ' << vert_2(D[n].coef * dta) / sqrt(n+1) << endl;
				if (vert_2(D[n].coef * dta) / sqrt(n+1) < eps) break;
				rhs = rhs - D[n].coef * dta;
			}
		}
		cout << "Dimension :" << 1 << endl;
		cout << "Condition :" << Cond_type << endl;
		cout << "Grid Number : " << n << endl;
		cout << "Restriction : " << (Restriction_type == Injection ? "Injection" : "Full Weighting") << endl;
		cout << "Interpolation : " << (Interpolation_type == Linear ? "Linear" : "Quadratic") << endl;
		cout << "Cycle : " << (Cycle_type == V_cycle ? "V_cycle" : "FMG") << endl;
		cout << "Iteration times : " << T1 << ", " << T2 << endl;
		cout << "Cycle times : " << iter << endl;
		Colvec<double> e = D[n].coef * sol - D[n].rhs;
		cout << "Residual Error in L_1 : " << vert_1(e) / (n+1) << endl;
		cout << "Residual Error in L_2 : " << vert_2(e) / sqrt(n+1) << endl;
		cout << "Residual Error in L_inf : " << vert_inf(e) << endl;
		if (is_test) {
			double h = 1.0 / n;
			double E1 = 0, E2 = 0, Einf = 0, C = 0;
			// Assume the average of solution be the arbitary constant.
			if (Cond_type == 3) {
				for (int i = 0; i <= n; ++ i) C += sol[i] - real(i*h);
				C /= (n+1);
			}
			for (int i = 0; i <= n; ++ i) {
				double ei = fabs(sol[i] - C - real(i*h));
				E1 += fabs(ei);
				E2 += ei * ei;
				Einf = max(Einf, ei);
			}
			E1 /= n+1, E2 /= n+1, E2 = sqrt(E2);
			cout << "Solution Error in L_1 : " << E1 << endl;
			cout << "Solution Error in L_2 : " << E2 << endl;
			cout << "Solution Error in L_inf : " << Einf << endl;
		}
	}
};

inline int id(int n, int i, int j) {
	return i * (n+1) + j;
}
Colvec<double> Injection_Restriction_2D(int n, const Colvec<double>& v) {
	Colvec<double> u((n/2+1)*(n/2+1));
	for (int i = 0; i <= n/2; ++ i)
		for (int j = 0; j <= n/2; ++ j)
			u[id(n/2, i, j)] = v[id(n, i*2, j*2)];
	return u;
}
Colvec<double> Full_Weighting_Restriction_2D(int n, const Colvec<double>& v) {
	Colvec<double> u((n/2+1)*(n/2+1));
	u[id(n/2,0,0)]		= v[id(n,0,0)];
	u[id(n/2,0,n/2)]	= v[id(n,0,n)];
	u[id(n/2,n/2,0)] 	= v[id(n,n,0)];
	u[id(n/2,n/2,n/2)]	= v[id(n,n,n)];
	for (int i = 1; i < n/2; ++ i) {
		u[id(n/2, i, 0)] = (
			  v[id(n, i*2, 0)]
			+ v[id(n, i*2-1, 0)] * 0.5
			+ v[id(n, i*2+1, 0)] * 0.5
		) / 2;
		u[id(n/2, i, n/2)] = (
			  v[id(n, i*2, n)]
			+ v[id(n, i*2-1, n)] * 0.5
			+ v[id(n, i*2+1, n)] * 0.5
		) / 2;
	}
	for (int j = 1; j < n/2; ++ j) {
		u[id(n/2, 0, j)] = (
			  v[id(n, 0, j*2)]
			+ v[id(n, 0, j*2-1)] * 0.5
			+ v[id(n, 0, j*2+1)] * 0.5
		) / 2;
		u[id(n/2, n/2, j)] = (
			  v[id(n, n, j*2)]
			+ v[id(n, n, j*2-1)] * 0.5
			+ v[id(n, n, j*2+1)] * 0.5
		) / 2;
	}
	for (int i = 1; i < n/2; ++ i)
		for (int j = 1; j < n/2; ++ j) {
			u[id(n/2, i, j)] = (
				  v[id(n, i*2, j*2)]
				+ v[id(n, i*2-1, j*2)] * 0.5
				+ v[id(n, i*2+1, j*2)] * 0.5
				+ v[id(n, i*2, j*2-1)] * 0.5
				+ v[id(n, i*2, j*2+1)] * 0.5
				+ v[id(n, i*2-1, j*2-1)] * 0.25
				+ v[id(n, i*2-1, j*2+1)] * 0.25
				+ v[id(n, i*2+1, j*2-1)] * 0.25
				+ v[id(n, i*2+1, j*2+1)] * 0.25
			) / 4;
		}
	return u;
}
Colvec<double> Linear_Interpolation_2D(int n, const Colvec<double>& v) {
	Colvec<double> u((n*2+1)*(n*2+1));
	for (int i = 0; i <= n; ++ i)
		for (int j = 0; j <= n; ++ j) {
			u[id(n*2, i*2, j*2)] = v[id(n, i, j)];
			if (i < n) u[id(n*2, i*2+1, j*2)] = (v[id(n, i, j)] + v[id(n, i+1, j)]) / 2;
			if (j < n) u[id(n*2, i*2, j*2+1)] = (v[id(n, i, j)] + v[id(n, i, j+1)]) / 2;
			if (i < n && j < n) u[id(n*2, i*2+1, j*2+1)] = (v[id(n, i, j)] + v[id(n, i+1, j)] + v[id(n, i, j+1)] + v[id(n, i+1, j+1)]) / 4;
		}
	return u;
}
Colvec<double> Quadratic_Interpolation_2D(int n, const Colvec<double>& v) {
	Colvec<double> u((n*2+1)*(n*2+1));
	for (int i = 0; i <= n; ++ i)
		for (int j = 0; j <= n; ++ j) {
			int i11, i12, i13, i21, i22, i23, i31, i32, i33;
			u[id(n*2, i*2, j*2)] = v[id(n, i, j)];
			if (i < n) {
				if (i < n/2)	i11 = id(n, i, j), i21 = id(n, i+1, j), i31 = id(n, i+2, j);
				else 			i11 = id(n, i+1, j), i21 = id(n, i, j), i31 = id(n, i-1, j);
				u[id(n*2, i*2+1, j*2)] = (v[i11] * 3 + v[i21] * 6 - v[i31]) / 8;
			}
			if (j < n) {
				if (j < n/2)	i11 = id(n, i, j), i12 = id(n, i, j+1), i13 = id(n, i, j+2);
				else			i11 = id(n, i, j), i12 = id(n, i, j-1), i13 = id(n, i, j-2);
				u[id(n*2, i*2, j*2+1)] = (v[i11] * 3 + v[i12] * 6 - v[i13]) / 8;
			}
			if (i < n && j < n) {
				if (i < n/2 && j < n/2) {
					i11 = id(n, i  , j), i12 = id(n, i  , j+1), i13 = id(n, i  , j+2);
					i21 = id(n, i+1, j), i22 = id(n, i+1, j+1), i23 = id(n, i+1, j+2);
					i31 = id(n, i+2, j), i32 = id(n, i+2, j+1), i33 = id(n, i+2, j+2);
				}
				else if (i < n/2 && j >= n/2) {
					i11 = id(n, i  , j+1), i12 = id(n, i  , j), i13 = id(n, i  , j-1);
					i21 = id(n, i+1, j+1), i22 = id(n, i+1, j), i23 = id(n, i+1, j-1);
					i31 = id(n, i+2, j+1), i32 = id(n, i+2, j), i33 = id(n, i+2, j-1);
				}
				else if (i >= n/2 && j < n/2) {
					i11 = id(n, i+1, j), i12 = id(n, i+1, j+1), i13 = id(n, i+1, j+2);
					i21 = id(n, i  , j), i22 = id(n, i  , j+1), i23 = id(n, i  , j+2);
					i31 = id(n, i-1, j), i32 = id(n, i-1, j+1), i33 = id(n, i-1, j+2);
				}
				else {
					i11 = id(n, i+1, j+1), i12 = id(n, i+1, j), i13 = id(n, i+1, j-1);
					i21 = id(n, i  , j+1), i22 = id(n, i  , j), i23 = id(n, i  , j-1);
					i31 = id(n, i-1, j+1), i32 = id(n, i-1, j), i33 = id(n, i-1, j-1);
				}
				u[id(n*2, i*2+1, j*2+1)] = (
					  v[i11] *  9 + v[i12] * 18 + v[i13] * -3
					+ v[i21] * 18 + v[i22] * 36 + v[i23] * -6
					+ v[i31] * -3 + v[i32] * -6 + v[i33] *  1
				) / 64;
			}
		}
	return u;
}


template <Cond_t Cond_type>
class Multigrid <2, Cond_type> {
private:
	double w;									// w is the release coefficient.
	Restriction_method Restriction_type;
	Interpolation_method Interpolation_type;
	Colvec<double> (*Restriction) (int, const Colvec<double>&);		// Restriction: Full Weighting and Injection.
	Colvec<double> (*Interpolation) (int ,const Colvec<double>&);	// Interpolation: Linear and Quadratic.
	Cycle_method Cycle_type;					// Cycle: V-cycle and FMG.
	map<int, Discretor<2, Cond_type>> D;			// Discretors for different grids.
	const Function_2D<double>& f;					// The rhs function.
	const Function_2D<double>& g;					// The boundary function.
	Colvec<double> sol;							// The solution vector.
public:
	Multigrid(const Function_2D<double>& f, const Function_2D<double>& g, Restriction_method Restriction_type, Interpolation_method Interpolation_type, Cycle_method Cycle_type):
		w(2.0/3), f(f), g(g), Restriction_type(Restriction_type), Interpolation_type(Interpolation_type), Cycle_type(Cycle_type) {
			if (Restriction_type == Injection) Restriction = Injection_Restriction_2D;
			else if (Restriction_type == Full_Weighting) Restriction = Full_Weighting_Restriction_2D;
			if (Interpolation_type == Linear) Interpolation = Linear_Interpolation_2D;
			else if (Interpolation_type == Quadratic) Interpolation = Quadratic_Interpolation_2D;
		}
private:
	Colvec<double> Jacobi(int n, Sparsed_Matrix<double>& A, Colvec<double>& b, const Colvec<double>& v0, int T) {
		int r = (n+1)*(n+1);
		Colvec<double> v(v0);
		Sparsed_Matrix<double> trns(r, r);
		Colvec<double> c(r);
		for (int i = 0; i < r; ++ i) {
			double dii = A[i][i];
			for (auto & [j, x] : A[i])
				if (j != i) trns[i][j] = -x / dii;
			c[i] = b[i] / dii;
		}
		Colvec<double> u(r);
		for (int i = 0; i < T; ++ i) {
			u = trns * v + c;
			v = w * u + (1-w) * v;
		}
		return v;
	}

	// v: Initial guess.
	// A: Discrete matrix for grid number n.
	// b: Discrete rhs for grid number n. 
	Colvec<double> VC(int n, Sparsed_Matrix<double>& A, Colvec<double>& b, Colvec<double>& v0, int T1, int T2) {
		int r = (n+1)*(n+1);
		Colvec<double> v = Jacobi(n, A, b, v0, T1);
		if (n <= 2) return Jacobi(n, A, b, v, T2);
		Colvec<double> c = Restriction(n, b - A * v);
		if (!D.count(n/2)) D.insert({n/2, Discretor<2, Cond_type>(n/2, f, g)});
		Colvec<double> zero((n/2+1)*(n/2+1));
		Colvec<double> v1 = VC(n/2, D[n/2].coef, c, zero, T1, T2);
		v = v + Interpolation(n/2, v1);
		return Jacobi(n, A, b, v, T2);
	}
	Colvec<double> FMGC(int n, Sparsed_Matrix<double>& A, Colvec<double>& b, int T1, int T2) {
		Colvec<double> zero((n+1)*(n+1));
		if (n <= 2) return VC(n, A, b, zero, T1, T2);
		Colvec<double> c = Restriction(n, b);
		if (!D.count(n/2)) D.insert({n/2, Discretor<2, Cond_type>(n/2, f, g)});
		Colvec<double> v = Interpolation(n/2, FMGC(n/2, D[n/2].coef, c, T1, T2));
		return VC(n, A, b, v, T1, T2);
	}
public:
	// n, the finest grid, must be a power of 2.
	// T1 and T2 are the Two-Grid iteration times specified by the user.
	// eps is the 
	void Solve(int n, int T1 = 5, int T2 = 5, double eps = 1e-8, bool is_test = 0, const Function_2D<double>& real = _0_2D<double>()) {
		if (!D.count(n)) D.insert({n, Discretor<2, Cond_type>(n, f, g)});
		Colvec<double> zero((n+1)*(n+1));
		Colvec<double> rhs = D[n].rhs;
		int r = (n+1)*(n+1);
		sol = zero;
		int iter = 0;
		if (Cycle_type == V_cycle) {
			while (1) {
				++ iter;
				Colvec<double> dta = VC(n, D[n].coef, rhs, zero, T1, T2);
				sol = sol + dta;
				cout << iter << ' ' << vert_2(D[n].coef * dta) / (n+1) << endl;
				if (vert_2(D[n].coef * dta) / (n+1) < eps) break;
				rhs = rhs - D[n].coef * dta;
			}
		}
		else {
			while (1) {
				++ iter;
				Colvec<double> dta = FMGC(n, D[n].coef, rhs, T1, T2);
				sol = sol + dta;
				cout << iter << ' ' << vert_2(D[n].coef * dta) / (n+1) << endl;
				if (vert_2(D[n].coef * dta) / (n+1) < eps) break;
				rhs = rhs - D[n].coef * dta;
			}
		}
		cout << "Dimension :" << 1 << endl;
		cout << "Condition :" << Cond_type << endl;
		cout << "Grid Number : " << n << endl;
		cout << "Restriction : " << (Restriction_type == Injection ? "Injection" : "Full Weighting") << endl;
		cout << "Interpolation : " << (Interpolation_type == Linear ? "Linear" : "Quadratic") << endl;
		cout << "Cycle : " << (Cycle_type == V_cycle ? "V_cycle" : "FMG") << endl;
		cout << "Iteration times : " << T1 << ", " << T2 << endl;
		cout << "Cycle times : " << iter << endl;
		Colvec<double> e = D[n].coef * sol - D[n].rhs;
		cout << "Residual Error in L_1 : " << vert_1(e) / r << endl;
		cout << "Residual Error in L_2 : " << vert_2(e) / sqrt(r) << endl;
		cout << "Residual Error in L_inf : " << vert_inf(e) << endl;
		if (is_test) {
			double h = 1.0 / n;
			double E1 = 0, E2 = 0, Einf = 0, C = 0;
			if (Cond_type == 15) {
				for (int i = 0; i <= n; ++ i)
					for (int j = 0; j <= n; ++ j)
						C += sol[id(n,i,j)] - real(i*h, j*h);
				C /= r;
			}
			for (int i = 0; i <= n; ++ i)
				for (int j = 0; j <= n; ++ j) if (i!=0&&i!=n || j!=0&&j!=n){
					double eij = fabs(sol[id(n,i,j)] - C - real(i*h, j*h));
					E1 += fabs(eij);
					E2 += eij * eij;
					Einf = max(Einf, eij);
					if (fabs(eij) > 0.01) cout << "gg " << i << ' ' << j << endl;
				}
			E1 /= r, E2 /= r, E2 = sqrt(E2);
			cout << "Solution Error in L_1 : " << E1 << endl;
			cout << "Solution Error in L_2 : " << E2 << endl;
			cout << "Solution Error in L_inf : " << Einf << endl;
		}
	}
};