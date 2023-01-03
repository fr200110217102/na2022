#include <bits/stdc++.h>
#include "../HW1/function.h"
#include "Thomas.h"
using namespace std;

template <class type>
class LinearSpline {
private:
	vector <type> x, f;	// 结点，函数值
public:
	LinearSpline(const vector<type>& x, const vector<type>& f) : x(x), f(f) {}
	type GetValue(const type& _x) {
		int i = upper_bound(x.begin(), x.end(), _x) - x.begin() - 1;
		type t = _x - x[i], c0 = f[i], c1 = (f[i+1] - f[i]) / (x[i+1] - x[i]);
		return c0 + c1 * t;
	}
};

template <class type>
LinearSpline<type> LinearSplineInterpolation(const Function <type>& f, const type& l, const type& r, const int& n) {
	vector <type> x(n+1), y(n+1);
	for(int i = 0; i <= n; ++ i) x[i] = l + (r - l) / n * i, y[i] = f(x[i]);
	return LinearSpline<type>(x, y);
}

template <class type>
class CubicSpline {
private:
	vector <type> x, f, m, M;	// 结点, 函数值, 一阶导, 二阶导
public:
	CubicSpline(const vector<type>& x, const vector<type>& f, const string& mode = "Natural", const type& m0 = 0, const type& mn = 0) : x(x), f(f) {
		int n = x.size() - 1;
		m.resize(n+1);
		vector <type> l(n+1), u(n+1);	// lambda, mu
		for (int i = 1; i <= n-1; ++ i) {
			l[i] = (x[i] - x[i-1]) / (x[i+1] - x[i-1]);
			u[i] = (x[i+1] - x[i]) / (x[i+1] - x[i-1]);
		}

		vector <type> dq1(n+1), dq2(n+2);	// 一阶差商, 二阶差商
		vector <type> a(n+1), b(n), c(n), y(n+1);

		if (mode == "Natural" && (m0 != 0 || mn != 0)) throw "Invalid Parameter!";
		if (mode == "Complete" || mode == "Natural") {
			// 计算差商
			for (int i = 1; i <= n; ++ i)
				dq1[i] = (f[i] - f[i-1]) / (x[i] - x[i-1]);
			dq2[1] = (dq1[1] - m0) / (x[1] - x[0]);
			for (int i = 2; i <= n; ++ i)
				dq2[i] = (dq1[i] - dq1[i-1]) / (x[i] - x[i-2]);
			dq2[n+1] = (mn - dq1[n]) / (x[n] - x[n-1]);
			// 构造三对角方程组
			for (int i = 1; i <= n-1; ++ i) {
				a[i] = 2;
				b[i-1] = u[i];
				c[i] = l[i];
				y[i] = 6 * dq2[i+1];
			}
			a[0] = 2, c[0] = 1, y[0] = 6 * dq2[1];
			a[n] = 2, b[n-1] = 1, y[n] = 6 * dq2[n+1];
			// 调用 Thomas 算法解出结点处的二阶导 M
			M = Thomas(a, b, c, y);
			// 根据 M 的值计算 m
			m[0] = m0, m[n] = mn;
			for (int i = 1; i <= n-1; ++ i)
				m[i] = dq1[i+1] - (2 * M[i] + M[i+1]) * (x[i+1] - x[i]) / 6;
		}
		else if (mode == "Specified_Second_Derivatives") {
			// 计算差商
			for (int i = 1; i <= n; ++ i)
				dq1[i] = (f[i] - f[i-1]) / (x[i] - x[i-1]);
			for (int i = 2; i <= n; ++ i)
				dq2[i] = (dq1[i] - dq1[i-1]) / (x[i] - x[i-2]);
			// 构造三对角方程组
			for (int i = 1; i <= n-1; ++ i) {
				a[i] = 2;
				b[i-1] = u[i];
				c[i] = l[i];
				y[i] = 6 * dq2[i+1];
			}
			a[0] = 1, y[0] = m0;
			a[n] = 1, y[n] = mn;
			// 调用 Thomas 算法解出结点处的二阶导 M
			M = Thomas(a, b, c, y);
			// 根据 M 的值计算 m
			for (int i = 0; i <= n-1; ++ i)
				m[i] = dq1[i+1] - (2 * M[i] + M[i+1]) * (x[i+1] - x[i]) / 6;
			m[n] = dq1[n] - (2 * M[n] + M[n-1]) * (x[n-1] - x[n]) / 6;
		}
	}
	type GetValue(const type& _x) const {
		// 求出 x 所在结点区间，根据该区间上的解析式计算样条函数的值
		int i = upper_bound(x.begin(), x.end(), _x) - x.begin() - 1;
		type t = _x - x[i], c0 = f[i], c1 = m[i], c2 = M[i] / 2, c3 = (M[i+1] - M[i]) / (x[i+1] - x[i]) / 6;
		return c0 + t * (c1 + t * (c2 + c3 * t));
	}
};

template <class type>
CubicSpline<type> CubicSplineInterpolation(const Function <type>& f, const type& l, const type& r, const int& n, const string& mode = "Natural") {
	vector <type> x(n+1), y(n+1);
	for(int i = 0; i <= n; ++ i) x[i] = l + (r - l) / n * i, y[i] = f(x[i]);
	type m0, mn;
	if (mode == "Natural") m0 = mn = 0;
	else if (mode == "Complete") m0 = f.d(l), mn = f.d(r);
	else if (mode == "Specified_Second_Derivatives") m0 = f.d(l, 2), mn = f.d(r, 2);
	return CubicSpline<type>(x, y, mode, m0, mn);
}