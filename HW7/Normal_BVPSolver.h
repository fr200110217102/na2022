#include <bits/stdc++.h>
#include "function.h"
#include "../Matrix.h"
#include "BVPSolver.h"
using namespace std;

// 输入是 n,f,g。
// 要求 f 的定义域是单位矩形，g 的定义域是矩形的四条边。f 和 g 都二阶可微。
// n 是分辨率，即网格每条边的格点数。
// 输出是解 sol，即所有格点处的函数值。
class Normal_BVPSolver : public BVPSolver{
private:
	int n, d;
	double h;				// 网格宽度
	Function_2D<double> &f, &g;
	Matrix<double> coef;	// 矩阵
	Colvec<double> rhs;		// 右端项
	Colvec<double> sol;		// 解
	map<pnt, int> id;
	vector<pnt> ps;
	string cond;
	void ins(pnt p) {
		if(!id.count(p)) id[p] = d++, ps.push_back(p);
	}
	void ID_Generator() {
		for (int i = 0; i <= n+1; ++ i)
			for (int j = 0; j <= n+1; ++ j)
				if (i!=0&&i!=n+1 || j!=0&&j!=n+1)
					ins(pnt(i*h,j*h));
	}
	void Laplace_Normal_Discretor(pnt p0, pnt p1, pnt p2, pnt p3, pnt p4) {
		int i0 = id[p0], i1 = id[p1], i2 = id[p2], i3 = id[p3], i4 = id[p4];
		coef[i0][i0] += 4 / (h*h);
		coef[i0][i1] -= 1 / (h*h);
		coef[i0][i2] -= 1 / (h*h);
		coef[i0][i3] -= 1 / (h*h);
		coef[i0][i4] -= 1 / (h*h);
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
	void Matrix_Generator() {
		coef = Matrix<double> (d, d);
		rhs = Colvec<double> (d);
		for (int i = 1; i <= n; ++ i)
			for (int j = 1; j <= n; ++ j) 
				Laplace_Normal_Discretor(pnt(i*h,j*h), pnt((i-1)*h,j*h), pnt((i+1)*h,j*h), pnt(i*h,(j-1)*h), pnt(i*h,(j+1)*h));
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
		
		if (cond == "NNNN") coef[0][0] += 1;
	}
public:
	Normal_BVPSolver(int n, Function_2D<double>& f, Function_2D<double>& g, const string& s) :
		n(n), d(0), h(1.0/(n+1)), f(f), g(g), cond(s) {}
	void Solve() {
		ID_Generator();
		Matrix_Generator();
		sol = Gauss_Improved_Solve(coef, rhs);
	}
	void Summary(Function_2D<double>& u, bool detail = 0) {
		cout << "Domain : Normal" << endl;
		cout << "Condition : " << cond << endl;
		cout << "n = " << n << ", h = " << h << endl;
		if (detail) cout << "Values :" << endl;
		double C = 0;
		if (cond == "NNNN") {
			for (int i = 0; i < d; ++ i) C += sol[i] - u(ps[i].x, ps[i].y);
			C /= d;
		}
		double res1 = 0, res2 = 0, resm = 0;
		int cnt = 0;
		for (int i = 0; i <= n+1; ++ i)
			for (int j = 0; j <= n+1; ++ j) if (i!=0&&i!=n+1 || j!=0&&j!=n+1) {
				++cnt;
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