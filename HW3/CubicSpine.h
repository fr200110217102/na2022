#include <bits/stdc++.h>
#include "../HW1/function.h"
#include "../Matrix.h"
using namespace std;

template <class type>
class CubicSpine {
private:
	int n;
	vector <type> x, f, m, M;	// 结点, 函数值, 一阶导, 二阶导
public:
	CubicSpine(const vector<type>& x, const vector<type>& f, const type& m0, const type& mn) : x(x), f(f) {
		n = x.size() - 1;
		m.resize(n+1);
		M.resize(n+1);
		vector <type> l(n+1), u(n+1);	// lambda, mu
		for (int i = 1; i <= n-1; ++ i) {
			l[i] = (x[i] - x[i-1]) / (x[i+1] - x[i-1]);
			u[i] = (x[i+1] - x[i]) / (x[i+1] - x[i-1]);
		}
		vector <type> dq1(n+1), dq2(n+2);	// 一阶差商, 二阶差商
		Matrix <type> A(n+1, n+1);			// 系数矩阵
		Colvec <type> b(n+1), sol(n+1);

		// 解 m
		for (int i = 1; i <= n; ++ i)
			dq1[i] = (f[i] - f[i-1]) / (x[i] - x[i-1]);
		A[0][0] = 1, A[n][n] = 1, b[0] = m0, b[n] = mn;
		for (int i = 1; i <= n-1; ++ i) {
			A[i][i] = 2;
			A[i][i-1] = l[i];
			A[i][i+1] = u[i];
			b[i] = 3 * (u[i] * dq1[i+1] + l[i] * dq1[i]);
		}
		sol = Gauss_Improved_Solve(A, b);
		for (int i = 0; i <= n; ++ i) m[i] = sol[i];

		// 解 M
		dq2[1] = (dq1[1] - m0) / (x[1] - x[0]);
		for (int i = 2; i <= n; ++ i)
			dq2[i] = (dq1[i] - dq1[i-1]) / (x[i] - x[i-2]);
		dq2[n+1] = (mn - dq1[n]) / (x[n] - x[n-1]);
		for (int i = 1; i <= n-1; ++ i) {
			A[i][i] = 2;
			A[i][i-1] = u[i];
			A[i][i+1] = l[i];
			b[i] = 6 * dq2[i+1];
		}
		A[0][0] = 2, A[0][1] = 1, b[0] = 6 * dq2[1];
		A[n][n] = 2, A[n][n-1] = 1, b[n] = 6 * dq2[n+1];
		sol = Gauss_Improved_Solve(A, b);
		for (int i = 0; i <= n; ++ i) M[i] = sol[i];
	}
	type GetValue(const type& _x) const {
		int i = upper_bound(x.begin(), x.end(), _x) - x.begin() - 1;
		type t = _x - x[i], c0 = f[i], c1 = m[i], c2 = M[i] / 2, c3 = (M[i+1] - M[i]) / (x[i+1] - x[i]) / 6;
		return c0 + t * (c1 + t * (c2 + c3 * t));
	}
};

template <class type>
CubicSpine<type> CubicSpineInterpolation(const Function <type>& f, const type& l, const type& r, const int& n) {
	vector <type> x(n+1), y(n+1);
	for(int i = 0; i <= n; ++ i) x[i] = l + (r - l) / n * i, y[i] = f(x[i]);
	type m0 = f.d(l), mn = f.d(r);
	return CubicSpine<type>(x, y, m0, mn);
}