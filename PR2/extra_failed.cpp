#include <bits/stdc++.h>
#include "function.h"
#include "Matrix.h"
#include "Sparsed_Matrix.h"
#include "multigrid.h"
using namespace std;

const double pi = acos(-1);

inline int sgn(double x) {
	if (fabs(x) < 1e-12) return 0;
	return x>0 ? 1 : -1;
}

struct vec {
	double x,y;
	vec() {}
	vec(double x,double y):x(x),y(y){}
	vec operator + (const vec& p) const {
		return vec(x+p.x, y+p.y);
	}
	vec operator - (const vec& p) const {
		return vec(x-p.x, y-p.y);
	}
	double norm() const {
		return sqrt(x*x + y*y);
	}
	bool operator < (const vec& p) const {
		return sgn(x-p.x) == -1 || (sgn(x-p.x) == 0 && sgn(y-p.y) == -1);
	}
	bool operator == (const vec& p) const {
		return sgn(x-p.x) == 0 && sgn(y-p.y) == 0;
	}
};
typedef vec pnt;

inline double dis(const pnt& p, const pnt& q) {
	return (p-q).norm();
}

typedef unsigned int Cond_t;
template <Cond_t Cond_type>
class Irnormal_Discretor {
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
	void Irnormal_Dirichlet_Discretor_2D(pnt p0, int i0, int i1, int i2) {
		pnt p1 = pnt(i1/(n+1)*h, i1%(n+1)*h);
		pnt p2 = pnt(i2/(n+1)*h, i2%(n+1)*h);
		double c1 = (p0.y-p1.y)/h, c2 = 1-c1;
		coef[i0][i1] = c2;
		coef[i0][i2] = c1;
		rhs[i0] = g(p0.x, p0.y);
	}
	void Irnormal_Neumann_Discretor_2D(pnt p0, int i0, const vector<int>& id) {
		vector<pnt> p;
		for (int i:id) p.push_back(pnt(i/(n+1)*h, i%(n+1)*h));
		Matrix<double> w(6, 6);
		Colvec<double> b(6);
		for (int i = 0; i < 6; ++ i) {
			w[0][i] = 1;
			w[1][i] = p[i].x - p0.x;
			w[2][i] = p[i].y - p0.y;
			w[3][i] = (p[i].x - p0.x) * (p[i].x - p0.x) / 2;
			w[4][i] = (p[i].y - p0.y) * (p[i].y - p0.y) / 2;
			w[5][i] = (p[i].x - p0.x) * (p[i].y - p0.y);
		}
		double dx = pi/16*cos(pi*p0.x), dy = -1;
		double dr = sqrt(dx*dx + dy*dy);
		b[1] = dx/dr, b[2] = dy/dr;
		Colvec<double> sol = Gauss_Improved_Solve(w, b);
		for (int i = 0; i < 6; ++ i) coef[i0][id[i]] = sol[i];
		rhs[i0] = g(p0.x, p0.y);
	}
public:
	Irnormal_Discretor() : n(0), f(_0_2D<double>()), g(_0_2D<double>()), coef(0, 0), rhs(0) {}
	Irnormal_Discretor(int n, const Function_2D<double>& f, const Function_2D<double>& g)
		:n(n), h(1.0/n), f(f), g(g), coef((n+1)*(n+1), (n+1)*(n+1)), rhs((n+1)*(n+1)) {
		h = 1.0/n;
		for (int i = 1; i < n; ++ i) {
			if (Cond_type & 1) {
				pnt p0(i*h, sin(pi*i*h)/16);
				int j = int(p0.y/h);
				if (i < n/2)	Irnormal_Neumann_Discretor_2D(p0, i, {id(i,j), id(i,j+1), id(i,j+2), id(i+1,j), id(i+1,j+1), id(i+2,j)});
				else			Irnormal_Neumann_Discretor_2D(p0, i, {id(i,j), id(i,j+1), id(i,j+2), id(i+1,j), id(i+1,j+1), id(i+2,j)});
				//Normal_Neumann_Discretor_2D(id(i,0), id(i,1), id(i,2));
			}
			else {
				pnt p0(i*h, sin(pi*i*h)/16);
				int j = int(p0.y/h);
				Irnormal_Dirichlet_Discretor_2D(p0, id(i,0), id(i,j), id(i,j+1));
			}
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
	}
};

template<Cond_t Cond_type>
class Irnormal_Multigrid {
private:
	double w;									// w is the release coefficient.
	Colvec<double> (*Restriction) (int, const Colvec<double>&);		// Restriction: Full Weighting and Injection.
	Colvec<double> (*Interpolation) (int ,const Colvec<double>&);	// Interpolation: Linear and Quadratic.
	Cycle_method Cycle_type;					// Cycle: V-cycle and FMG.
	map<int, Irnormal_Discretor<Cond_type>> D;	// Discretors for different grids.
	const Function_2D<double>& f;				// The rhs function.
	const Function_2D<double>& g;				// The boundary function.
	Colvec<double> sol;							// The solution vector.
public:
	Irnormal_Multigrid(const Function_2D<double>& f, const Function_2D<double>& g): w(2.0/3), f(f), g(g) {}
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
		if (!D.count(n/2)) D.insert({n/2, Irnormal_Discretor<Cond_type>(n/2, f, g)});
		Colvec<double> zero((n/2+1)*(n/2+1));
		Colvec<double> v1 = VC(n/2, D[n/2].coef, c, zero, T1, T2);
		v = v + Interpolation(n/2, v1);
		return Jacobi(n, A, b, v, T2);
	}
	Colvec<double> FMGC(int n, Sparsed_Matrix<double>& A, Colvec<double>& b, int T1, int T2) {
		Colvec<double> zero((n+1)*(n+1));
		if (n <= 4) return VC(n, A, b, zero, T1, T2);
		Colvec<double> c = Restriction(n, b);
		if (!D.count(n/2)) D.insert({n/2, Irnormal_Discretor<Cond_type>(n/2, f, g)});
		Colvec<double> v = Interpolation(n/2, FMGC(n/2, D[n/2].coef, c, T1, T2));
		return VC(n, A, b, v, T1, T2);
	}
public:
	// n, the finest grid, must be a power of 2.
	// T1 and T2 are the Two-Grid iteration times specified by the user.
	void Solve(int n, Restriction_method Restriction_type, Interpolation_method Interpolation_type, Cycle_method Cycle_type, int T1 = 5, int T2 = 5, double eps = 1e-8, bool is_test = 0, const Function_2D<double>& real = _0_2D<double>()) {
		if (Restriction_type == Injection) Restriction = Injection_Restriction_2D;
		else if (Restriction_type == Full_Weighting) Restriction = Full_Weighting_Restriction_2D;
		if (Interpolation_type == Linear) Interpolation = Linear_Interpolation_2D;
		else if (Interpolation_type == Quadratic) Interpolation = Quadratic_Interpolation_2D;
		if (!D.count(n)) D.insert({n, Irnormal_Discretor<Cond_type>(n, f, g)});
		Colvec<double> zero((n+1)*(n+1));
		Colvec<double> rhs = D[n].rhs;
		Colvec<double> sol = Gauss_Improved_Solve((Matrix<double>)D[n].coef, D[n].rhs);
		int r = (n+1)*(n+1);
/*		sol = zero;
		int iter = 0;
		double start_time = clock();
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
		double end_time = clock();*/
		cout << "Dimension :" << 2 << endl;
		cout << "Condition :" << Cond_type << endl;
		cout << "Grid Number : " << n << endl;
		cout << "Restriction : " << (Restriction_type == Injection ? "Injection" : "Full Weighting") << endl;
		cout << "Interpolation : " << (Interpolation_type == Linear ? "Linear" : "Quadratic") << endl;
		cout << "Cycle : " << (Cycle_type == V_cycle ? "V_cycle" : "FMG") << endl;
		cout << "Iteration times : " << T1 << ", " << T2 << endl;
//		cout << "Cycle times : " << iter << endl;
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
				}
			E1 /= r, E2 /= r, E2 = sqrt(E2);
			cout << "Solution Error in L_1 : " << E1 << endl;
			cout << "Solution Error in L_2 : " << E2 << endl;
			cout << "Solution Error in L_inf : " << Einf << endl;
		}
//		cout << "Run time : " << (end_time - start_time) / CLOCKS_PER_SEC << endl;
	}
};

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
		if (x == 0) return Cond_type & 2 ? u.partial(0, y, 1, 0) : u(0, y) ;
		if (y == 1) return Cond_type & 4 ? -u.partial(x, 1, 0, 1) : u(x, 1);
		if (x == 1) return Cond_type & 8 ? -u.partial(1, y, 1, 0) : u(1, y);
		if (fabs(y - sin(pi * x) / 16) < eps) {
			if (Cond_type & 1) {
				double dx = pi*cos(pi*x)/16, dy = -1;
				double dr = sqrt(dx*dx + dy*dy);
				return (u.partial(x, y, 1, 0) * dx + u.partial(x, y, 0, 1) * dy) / dr;
			}
			else return u(x, y);
		}
		cout << "?" << x<<' '<<y<<endl;
		throw 0;
	}
};

int main(int argc, char** argv) {
	int n = atoi(argv[1]);
	const Restriction_method i1 = string(argv[2]) == "I" ? Injection : Full_Weighting;
	const Interpolation_method i2 = string(argv[3]) == "L" ? Linear : Quadratic;
	const Cycle_method i3 = string(argv[4]) == "V" ? V_cycle : FMG;
	int T1 = atoi(argv[5]), T2 = atoi(argv[6]);
	auto Solver_D = Irnormal_Multigrid<0>(F(u), G<0>(u));
	Solver_D.Solve(n, i1, i2, i3, T1, T2, 1e-8, 1, u);
	auto Solver_N = Irnormal_Multigrid<15>(F(u), G<15>(u));
	Solver_N.Solve(n, i1, i2, i3, T1, T2, 1e-8, 1, u);
	auto Solver_M = Irnormal_Multigrid<3>(F(u), G<3>(u));
	Solver_M.Solve(n, i1, i2, i3, T1, T2, 1e-8, 1, u);
}