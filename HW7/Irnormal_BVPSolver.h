#include <bits/stdc++.h>
#include "function.h"
#include "geo_2D.h"
#include "../Matrix.h"
#include "BVPSolver.h"
using namespace std;

class Irnormal_BVPSolver : public BVPSolver{
private:
	int n, d;
	circle D;
	double h;
	Function_2D<double> &f, &g;
	Matrix<double> coef;
	Colvec<double> rhs;
	Colvec<double> sol;
	map<pnt, int> id;
	vector<pnt> ps;
	string cond;
	void ins(pnt p) {
		if(!id.count(p)) id[p] = d++, ps.push_back(p);
	}
	void ID_Generator() {
		for (int i = 1; i <= n; ++ i)
			for (int j = 1; j <= n; ++ j)
				if (pnt_to_circle(pnt(i*h, j*h), D) >= 0) 
					ins(pnt(i*h, j*h));
		for (int i = 1; i <= n ;++ i) {
			ins(pnt(0, i*h));
			ins(pnt(1, i*h));
			ins(pnt(i*h, 0));
			ins(pnt(i*h, 1));
		}
		for (int i = 1; i <= n; ++ i)
			if (sgn(i*h - (D.o.x-D.r)) * sgn(i*h - (D.o.x+D.r)) <= 0) {
				ins(pnt(i*h, D.Y1(i*h)));
				ins(pnt(i*h, D.Y2(i*h)));
			}
		for (int j = 1; j <= n; ++ j)
			if (sgn(j*h - (D.o.y-D.r)) * sgn(j*h - (D.o.y+D.r)) <= 0) {
				ins(pnt(D.X1(j*h), j*h));
				ins(pnt(D.X2(j*h), j*h));
			}
	}
	void Laplace_Irnormal_Discretor(pnt p0, pnt p1, pnt p2, pnt p3, pnt p4) {
		int i0 = id[p0], i1 = id[p1], i2 = id[p2], i3 = id[p3], i4 = id[p4];
		double lt = (p0.x - p1.x) / h, rt = (p2.x - p0.x) / h;
		double dt = (p0.y - p3.y) / h, ut = (p4.y - p0.y) / h;
		coef[i0][i0] += (lt+rt) / ((lt+rt)*lt*rt*h*h/2);
		coef[i0][i1] -= rt / ((lt+rt)*lt*rt*h*h/2);
		coef[i0][i2] -= lt / ((lt+rt)*lt*rt*h*h/2);
		coef[i0][i0] += (dt+ut) / ((dt+ut)*dt*ut*h*h/2);
		coef[i0][i3] -= ut / ((dt+ut)*dt*ut*h*h/2);
		coef[i0][i4] -= dt / ((dt+ut)*dt*ut*h*h/2);
		rhs[i0] = f(p0.x, p0.y);
	}
	void Dirichlet_Discretor(pnt p0) {
		int i0 = id[p0];
		coef[i0][i0] = 1;
		rhs[i0] = g(p0.x, p0.y);
	}
	void Neumann_Normal_Discretor(pnt p0, pnt p1, pnt p2) {
		int i0 = id[p0], i1 = id[p1], i2 = id[p2];
		coef[i0][i0] -= 1.5 / h;
		coef[i0][i1] += 2 / h;
		coef[i0][i2] -= 0.5 / h;
		rhs[i0] = g(p0.x, p0.y);
	}
	void Neumann_Irnormal_Discretor(pnt p0, pnt p1, pnt p2, pnt p3, pnt p4, pnt p5) {
		int i0 = id[p0], i1 = id[p1], i2 = id[p2], i3 = id[p3], i4 = id[p4], i5 = id[p5];
		Matrix<double> c(6, 6);
		Colvec<double> b(6);
		pnt p[6] = {p0, p1, p2, p3, p4, p5};
		for (int i = 0; i < 6; ++ i) {
			c[0][i] = 1;
			c[1][i] = p[i].x-p0.x;
			c[2][i] = p[i].y-p0.y;
			c[3][i] = (p[i].x-p0.x) * (p[i].x-p0.x) / 2;
			c[4][i] = (p[i].y-p0.y) * (p[i].y-p0.y) / 2;
			c[5][i] = (p[i].x-p0.x) * (p[i].y-p0.y);
		}
		double arg = D.arg(p0);
		b[1] = cos(arg);
		b[2] = sin(arg);
		Colvec<double> w = Gauss_Improved_Solve(c, b);
		coef[i0][i0] = w[0];
		coef[i0][i1] = w[1];
		coef[i0][i2] = w[2];
		coef[i0][i3] = w[3];
		coef[i0][i4] = w[4];
		coef[i0][i5] = w[5];
		rhs[i0] = g(p0.x, p0.y);
	}
	void Matrix_Generator() {
		coef = Matrix<double> (d, d);
		rhs = Colvec<double> (d);
		for (int i = 1; i <= n; ++ i)
			for (int j = 1; j <= n; ++ j) if (pnt_to_circle(pnt(i*h, j*h), D) > 0){
				pnt p(i*h, j*h);
				pnt lp((i-1)*h, j*h);
				pnt rp((i+1)*h, j*h);
				pnt dp(i*h, (j-1)*h);
				pnt up(i*h, (j+1)*h);
				if (pnt_to_circle(lp, D) < 0) lp.x = D.X2(j*h);
				if (pnt_to_circle(rp, D) < 0) rp.x = D.X1(j*h);
				if (pnt_to_circle(dp, D) < 0) dp.y = D.Y2(i*h);
				if (pnt_to_circle(up, D) < 0) up.y = D.Y1(i*h);
				Laplace_Irnormal_Discretor(p, lp, rp, dp, up);
			}
		
		for (int i = 1; i <= n; ++ i)
			if (cond[0] == 'D') Dirichlet_Discretor(pnt(i*h, 0));
			else Neumann_Normal_Discretor(pnt(i*h, 0), pnt(i*h, h), pnt(i*h, 2*h));

		for (int j = 1; j <= n; ++ j)
			if (cond[1] == 'D') Dirichlet_Discretor(pnt(0, j*h));
			else Neumann_Normal_Discretor(pnt(0, j*h), pnt(h, j*h), pnt(2*h, j*h));

		for (int i = 1; i <= n; ++ i)
			if (cond[2] == 'D') Dirichlet_Discretor(pnt(i*h, 1));
			else Neumann_Normal_Discretor(pnt(i*h, 1), pnt(i*h, 1-h), pnt(i*h, 1-2*h));
			
		for (int j = 1; j <= n; ++ j)
			if (cond[3] == 'D') Dirichlet_Discretor(pnt(1, j*h));
			else Neumann_Normal_Discretor(pnt(1, j*h), pnt(1-h, j*h), pnt(1-2*h, j*h));
		
		if (cond[4] == 'D') {
			for (int i = 1; i <= n; ++ i) 
				if (sgn(i*h - (D.o.x-D.r)) * sgn(i*h - (D.o.x+D.r)) == -1) {
					Dirichlet_Discretor(pnt(i*h, D.Y1(i*h)));
					Dirichlet_Discretor(pnt(i*h, D.Y2(i*h)));
				}
			for (int j = 1; j <= n; ++ j)
				if (sgn(j*h - (D.o.y-D.r)) * sgn(j*h - (D.o.y+D.r)) == -1) {
					Dirichlet_Discretor(pnt(D.X1(j*h), j*h));
					Dirichlet_Discretor(pnt(D.X2(j*h), j*h));
				}
		}
		else {
			for (int i = 1, j; i <= n; ++ i)
				if (sgn(i*h - (D.o.x-D.r)) * sgn(i*h - (D.o.x+D.r)) <= 0) {
					pnt p(i*h, D.Y1(i*h));
					int j = int(p.y/h-1e-12);
					int op = sgn(i*h-D.o.x);
					if(op == 0) op = 1;
					pnt p1(i*h, j*h);
					pnt p2(i*h, (j-1)*h);
					pnt p3((i+op)*h, j*h);
					pnt p4((i+op)*h, (j-1)*h);
					pnt p5((i+2*op)*h, j*h);
					Neumann_Irnormal_Discretor(p, p1, p2, p3, p4, p5);
				}
			for (int i = 1; i <= n; ++ i)
				if (sgn(i*h - (D.o.x-D.r)) * sgn(i*h - (D.o.x+D.r)) == -1) {
					pnt p(i*h, D.Y2(i*h));
					int j = int(p.y/h+1e-12) + 1;
					int op = sgn(i*h-D.o.x);
					if (op == 0) op = 1;
					pnt p1(i*h, j*h);
					pnt p2(i*h, (j+1)*h);
					pnt p3((i+op)*h, j*h);
					pnt p4((i+op)*h, (j+1)*h);
					pnt p5((i+2*op)*h, j*h);
					Neumann_Irnormal_Discretor(p, p1, p2, p3, p4, p5);
				}
			for (int j = 1; j <= n; ++ j)
				if (sgn(j*h - (D.o.y-D.r)) * sgn(j*h - (D.o.y+D.r)) <= 0) {
					pnt p(D.X1(j*h), j*h);
					int i = int(p.x/h-1e-12);
					int op = sgn(j*h-D.o.y);
					if (op == 0) op = 1;
					pnt p1(i*h, j*h);
					pnt p2((i-1)*h, j*h);
					pnt p3(i*h, (j+op)*h);
					pnt p4((i-1)*h, (j+op)*h);
					pnt p5(i*h, (j+2*op)*h);
					Neumann_Irnormal_Discretor(p, p1, p2, p3, p4, p5);
				}
			for (int j = 1; j <= n; ++ j)
				if (sgn(j*h - (D.o.y-D.r)) * sgn(j*h - (D.o.y+D.r)) == -1) {
					pnt p(D.X2(j*h), j*h);
					int i = int(p.x/h+1e-12) + 1;
					int op = sgn(j*h-D.o.y);
					if (op == 0) op = 1;
					pnt p1(i*h, j*h);
					pnt p2((i+1)*h, j*h);
					pnt p3(i*h, (j+op)*h);
					pnt p4((i+1)*h, (j+op)*h);
					pnt p5(i*h, (j+2*op)*h);
					Neumann_Irnormal_Discretor(p, p1, p2, p3, p4, p5);
				}
		}
		if (cond == "NNNNN") coef[0][0] += 1;
	}
public:
	Irnormal_BVPSolver(int n, Function_2D<double>& f, Function_2D<double>& g, double x0, double y0, double r, const string& s) :
		n(n), d(0), ps(0), h(1.0/(n+1)), f(f), g(g), D(x0, y0, r), cond(s) {
			if (r <= h) {
				cerr << "Too Coarse Grid!" << endl;
				throw 1;
			}
			if (x0 - r <= 2*h || x0 + r >= 1-2*h || y0 - r <= 2*h || y0 + r >= 1-2*h) {
				cerr << "Circle Outside the Square!" << endl;
				throw 1;
			}
		}
	void Solve() {
		ID_Generator();
		Matrix_Generator();
		sol = Gauss_Improved_Solve(coef, rhs);
	}
	void Summary(Function_2D<double>& u, bool detail = 0) {
		cout << "Domain : Irnormal" << endl;
		cout << "Condition : " << cond << endl;
		cout << "n = " << n << ", h = " << h << endl;
		cout << "Circle : " << "O = (" << D.o.x << ", " << D.o.y << "), r = " << D.r << endl;
		if (detail) cout << "Values :" << endl;
		double C = 0;
		if (cond == "NNNNN") {
			for (int i = 0; i < d; ++ i) C += sol[i] - u(ps[i].x, ps[i].y);
			C /= d;
		}
		double res1 = 0, res2 = 0, resm = 0;
		int cnt = 0;
		for (int i = 1; i <= n; ++ i)
			for (int j = 1; j <= n; ++ j) if (pnt_to_circle(pnt(i*h, j*h), D) > 0) {
				++ cnt;
				double uij = sol[id[pnt(i*h,j*h)]], uij_real = u(i*h, j*h);
				double eij = fabs(uij - C - uij_real);
				res1 += eij;
				res2 += eij*eij;
				resm = max(resm, eij);
				if (detail) cout << "(" << i*h << ", " << j*h << "), " << "Solution Value : " << uij << ", Real Value : " << uij_real << endl;
			}
		res1 /= cnt, res2 /= cnt, res2 = sqrt(res2);
		cout << "Solution Error In L_1 : " << res1 << endl;
		cout << "Solution Error In L_2 : " << res2 << endl;
		cout << "Solution Error In L_max : " << resm << endl;
	}
};