#include <bits/stdc++.h>
#include "../HW1/function.h"
#include "../Matrix.h"
using namespace std;

template <class type>
type B3(const int& i, const type& x) {
	if (x > i-1 && x <= i) return (x-i+1) * (x-i+1) * (x-i+1) / 6;
	if (x > i && x <= i+1) return 2.0 / 3 - (x-i+1) * (i+1-x) * (i+1-x) / 2;
	if (x > i+1 && x <= i+3) return B3(i, 2*i+2-x);
	return 0;
}

template <class type>
class CubicBSpline{
private:
	int n, L, R;
	vector <type> f, a;
public:
	CubicBSpline(const int& x0, const int& n, const vector<type>& f, const type& m0, const type& mn) : n(n), L(x0), f(f) {
		R = L + n;
		Matrix <type> M(n+1, n+1);
		Colvec <type> b(n+1);
		for (int i = 1; i < n; ++ i) {
			M[i][i] = 4;
			M[i][i-1] = 1;
			M[i][i+1] = 1;
			b[i] = 6 * f[i];
		}
		M[0][0] = 2, M[0][1] = 1, b[0] = 3 * f[0] + m0;
		M[n][n] = 2, M[n][n-1] = 1, b[n] = 3 * f[n] - mn;
		Colvec <type> sol = Cholesky_Improved_Solve(M, b);
		a.resize(n+3);
		for (int i = 0; i <= n; ++ i) a[i+1] = sol[i];
		a[0] = a[2] - 2 * m0;
		a[n+2] = a[n] + 2 * mn;
	}
	type GetValue(const type& _x) const {
		if (_x < L || _x > R) return nan("");
		int i = floor(_x);
		type res = 0;
		if(i-2 >= L-2 && i-2 <= R) res += a[i-2-L+2] * B3(i-2, _x);
		if(i-1 >= L-2 && i-1 <= R) res += a[i-1-L+2] * B3(i-1, _x);
		if(i   >= L-2 && i   <= R) res += a[i  -L+2] * B3(i  , _x);
		if(i+1 >= L-2 && i+1 <= R) res += a[i+1-L+2] * B3(i+1, _x);
		return res;
	}
};

template <class type>
CubicBSpline<type> CubicBSplineInterpolation(const Function <type>& f, const int& l, const int& n) {
	vector <type> y(n+1);
	for(int i = 0; i <= n; ++ i) y[i] = f(l+i);
	type m0 = f.d(l), mn = f.d(l+n);
	return CubicBSpline<type>(l, n, y, m0, mn);
}

template <class type>
type B2(const int& i, const type& x) {
	if (x > i-1 && x <= i) return (x-i+1) * (x-i+1) / 2;
	if (x > i && x <= i+1) return 0.75 - (x-i-0.5) * (x-i-0.5);
	if (x > i+1 && x <= i+2) return (i+2-x) * (i+2-x) / 2;
	return 0;
}

template <class type>
class QuadraticBSpline{
private:
	int n, L, R;
	vector <type> f, a;
public:
	QuadraticBSpline(const int& x0, const int& n, const vector<type>& f, const type& f0, const type& fn) : n(n), L(x0), f(f) {
		R = L + n;
		Matrix <type> M(n, n);
		Colvec <type> b(n);
		for (int i = 1; i < n-1; ++ i) {
			M[i][i] = 6;
			M[i][i-1] = 1;
			M[i][i+1] = 1;
			b[i] = 8 * f[i];
		}
		M[0][0] = 5, M[0][1] = 1, b[0] = 8 * f[0] - 2 * f0;
		M[n-1][n-1] = 5, M[n-1][n-2] = 1, b[n-1] = 8 * f[n-1] - 2 * fn;
		Colvec <type> sol = Cholesky_Improved_Solve(M, b);
		a.resize(n+2);
		for (int i = 0; i < n; ++ i) a[i+1] = sol[i];
		a[0] = 2 * f0 - a[1];
		a[n+1] = 2 * fn - a[n];
	}
	type GetValue(const type& _x) const {
		if (_x < L || _x > R) return nan("");
		int i = floor(_x);
		type res = 0;
		if(i-1 >= L-1 && i-1 <= R) res += a[i-1-L+1] * B2(i-1, _x);
		if(i   >= L-1 && i   <= R) res += a[i  -L+1] * B2(i  , _x);
		if(i+1 >= L-1 && i+1 <= R) res += a[i+1-L+1] * B2(i+1, _x);
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