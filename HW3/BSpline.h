#include <bits/stdc++.h>
#include "../HW1/function.h"
#include "Thomas.h"
using namespace std;

template <class type>
type B(const int& n, const int& i, const type& x) {
	if (n == 0) {
		if (i-1 < x && x <= i) return 1;
		else return 0;
	}
	// 根据 B 样条基的递归定义计算
	return (x-i+1) * B(n-1, i, x) / n + (i+n-x) * B(n-1, i+1, x) / n;
}

template <class type>
class LinearBSpline{
private:
	int L, R;
	vector<type> coef;
public:
	LinearBSpline(const int& x0, const int& n, const vector<type>& f) : L(x0), coef(f) {
		R = L + n;
	}
	type GetValue(const type& _x) const{
		if (_x < L || _x > R) throw "Out of Range!";
		int i = floor(_x);
		type res = 0;
		if(i   >= L && i   <= R) res += coef[i  -L] * B(1, i,   _x);
		if(i+1 >= L && i+1 <= R) res += coef[i+1-L] * B(1, i+1, _x);
		return res;
	}
};

template <class type>
LinearBSpline<type> LinearBSplineInterpolation(const Function <type>& f, const int& l, const int& n) {
	vector <type> y(n+1);
	for(int i = 0; i <= n; ++ i) y[i] = f(l+i);
	return LinearBSpline<type>(l, n, y);
}

template <class type>
class CubicBSpline{
private:
	int L, R;
	vector <type> coef;
public:
	CubicBSpline(const int& x0, const int& n, const vector<type>& f, const string& mode = "Natural", const type& m0 = 0, const type& mn = 0) : L(x0) {
		R = L + n;
		vector <type> a(n+1), b(n), c(n), y(n+1);
		// 构造三对角方程组的中间 n-1 行
		for (int i = 1; i < n; ++ i) {
			a[i] = 4;
			b[i-1] = 1;
			c[i] = 1;
			y[i] = 6 * f[i];
		}
		// 根据边界条件构造首尾两行
		if (mode == "Natural" || mode == "Complete") {
			a[0] = 2, c[0] = 1, y[0] = 3 * f[0] + m0;
			a[n] = 2, b[n-1] = 1, y[n] = 3 * f[n] - mn;
			// 解出中间 n+1 个系数
			vector <type> t = Thomas(a, b, c, y);
			coef.resize(n+3);
			for (int i = 0; i <= n; ++ i) coef[i+1] = t[i];
			// 计算首尾两个系数
			coef[0] = coef[2] - 2 * m0;
			coef[n+2] = coef[n] + 2 * mn;
		}
		else if (mode == "Specified_Second_Derivatives") {
			a[0] = 6, y[0] = 6 * f[0] - m0;
			a[n] = 6, y[n] = 6 * f[n] - mn;
			// 解出中间 n+1 个系数
			vector <type> t = Thomas(a, b, c, y);
			coef.resize(n+3);
			for (int i = 0; i <= n; ++ i) coef[i+1] = t[i];
			// 计算首尾两个系数
			coef[0] = m0 - coef[2] + 2 * coef[1];
			coef[n+2] = mn - coef[n] + 2 * coef[n-1];
		}
	}
	type GetValue(const type& _x) const {
		if (_x < L || _x > R) throw "Out of Range!";
		int i = floor(_x);
		type res = 0;
		if(i-2 >= L-2 && i-2 <= R) res += coef[i-2-L+2] * B(3, i-2, _x);
		if(i-1 >= L-2 && i-1 <= R) res += coef[i-1-L+2] * B(3, i-1, _x);
		if(i   >= L-2 && i   <= R) res += coef[i  -L+2] * B(3, i  , _x);
		if(i+1 >= L-2 && i+1 <= R) res += coef[i+1-L+2] * B(3, i+1, _x);
		return res;
	}
};

template <class type>
CubicBSpline<type> CubicBSplineInterpolation(const Function <type>& f, const int& l, const int& n, const string& mode = "Natural") {
	vector <type> y(n+1);
	for(int i = 0; i <= n; ++ i) y[i] = f(l+i);
	type m0, mn;
	if (mode == "Natural") m0 = mn = 0;
	else if (mode == "Complete") m0 = f.d(l), mn = f.d(l+n);
	else if (mode == "Specified_Second_Derivatives") m0 = f.d(l, 2), mn = f.d(l+n, 2);
	return CubicBSpline<type>(l, n, y, mode, m0, mn);
}

template <class type>
class QuadraticBSpline{
private:
	int L, R;
	vector <type> coef;
public:
	QuadraticBSpline(const int& x0, const int& n, const vector<type>& f, const type& f0, const type& fn) : L(x0) {
		R = L + n;
		vector <type> a(n), b(n-1), c(n-1), y(n);
		// 构造三对角方程组
		for (int i = 1; i < n-1; ++ i) {
			a[i] = 6;
			b[i-1] = 1;
			c[i] = 1;
			y[i] = 8 * f[i];
		}
		a[0] = 5, c[0] = 1, y[0] = 8 * f[0] - 2 * f0;
		a[n-1] = 5, b[n-2] = 1, y[n-1] = 8 * f[n-1] - 2 * fn;
		// 解出中间 n 个系数
		vector <type> t = Thomas(a, b, c, y);
		coef.resize(n+2);
		for (int i = 0; i < n; ++ i) coef[i+1] = t[i];
		// 计算首尾两个系数
		coef[0] = 2 * f0 - coef[1];
		coef[n+1] = 2 * fn - coef[n];
	}
	type GetValue(const type& _x) const {
		if (_x < L || _x > R) throw "Out of Range!";
		int i = floor(_x);
		type res = 0;
		if(i-1 >= L-1 && i-1 <= R) res += coef[i-1-L+1] * B(2, i-1, _x);
		if(i   >= L-1 && i   <= R) res += coef[i  -L+1] * B(2, i  , _x);
		if(i+1 >= L-1 && i+1 <= R) res += coef[i+1-L+1] * B(2, i+1, _x);
		return res;
	}
};

template <class type>
QuadraticBSpline<type> QuadraticBSplineInterpolation(const Function <type>& f, const int& l, const int& n) {
	vector <type> y(n);
	for(int i = 0; i < n; ++ i) y[i] = f(l+i+0.5);
	type f0 = f(l), fn = f(l+n);
	return QuadraticBSpline<type>(l, n, y, f0, fn);
}