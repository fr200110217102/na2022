#include <bits/stdc++.h>
#include "function.h"
#include "../Matrix.h"
#include "BVPSolver.h"
using namespace std;

// 输入是 n,f,g。
// 要求 f 的定义域是单位矩形，g 的定义域是矩形的四条边。f 和 g 都二阶可微。
// n 是分辨率，即网格每条边的格点数。
// 输出是解 sol，即所有格点处的函数值。
class Dirichlet_Normal_Solver : public BVPSolver{
private:
	int n, d;
	double h;				// 网格宽度
	Function_2D<double> &f, &g;
	Matrix<double> coef;	// 矩阵
	Colvec<double> rhs;		// 右端项
	Colvec<double> sol;		// 解
	void Matrix_Generator() {
		for (int i = 1; i <= n; ++ i)
			for (int j = 1; j <= n; ++ j) {
				coef[ID(i,j)][ID(i,j)] = 4 / (h*h);
				coef[ID(i,j)][ID(i-1,j)] = -1 / (h*h);
				coef[ID(i,j)][ID(i+1,j)] = -1 / (h*h);
				coef[ID(i,j)][ID(i,j-1)] = -1 / (h*h);
				coef[ID(i,j)][ID(i,j+1)] = -1 / (h*h);
				rhs[ID(i,j)] = f(i*h, j*h);
			}
		for (int j = 1; j <= n; ++ j) {
			coef[ID(0,j)][ID(0,j)] = 1;
			rhs[ID(0,j)] = g(0, j*h);
			coef[ID(n+1,j)][ID(n+1,j)] = 1;
			rhs[ID(n+1,j)] = g(1, j*h);
		}
		for (int i = 1; i <= n; ++ i) {
			coef[ID(i,0)][ID(i,0)] = 1;
			rhs[ID(i,0)] = g(i*h, 0);
			coef[ID(i,n+1)][ID(i,n+1)] = 1;
			rhs[ID(i,n+1)] = g(i*h, 1);
		}
		coef[ID(0,0)][ID(0,0)] = 1;
		coef[ID(0,n+1)][ID(0,n+1)] = 1;
		coef[ID(n+1,0)][ID(n+1,0)] = 1;
		coef[ID(n+1,n+1)][ID(n+1,n+1)] = 1;
	}
public:
	Dirichlet_Normal_Solver(int n, Function_2D<double>& f, Function_2D<double>& g) :
		n(n), d((n+2)*(n+2)), h(1.0/(n+1)), f(f), g(g), coef(d,d), rhs(d) {}
	// 对网格结点进行编号。
	int ID(int i, int j) {
		return i*(n+2) + j;
	}
	void Solve() {
		Matrix_Generator();
		sol = Gauss_Improved_Solve(coef, rhs);
		// sol = Gauss_Seidel(coef, rhs);
	}
	double getval(int i, int j) {
		return sol[ID(i,j)];
	}
};