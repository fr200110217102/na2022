#include <bits/stdc++.h>
#include "Spline.h"
#include "../HW1/function.h"
using namespace std;

const double PI = acos(-1), Q = sqrt(3), dt = 1e-8, eps = 1e-12;

class X : public Function<double> {
public:
	virtual double operator () (const double& t) const {
		return Q * sin(t);
	}
	virtual double d(const double& t) const {
		return Q * cos(t);
	}
} x;

class Y : public Function<double>{
public:
	virtual double operator () (const double& t) const {
		return 2 * (Q * cos(t) + sqrt(Q * fabs(sin(t)))) / 3;
	}
	virtual double d(const double& t) const {
		if (sin(t) > 0) return 2 * (-Q * sin(t) + sqrt(Q) * cos(t) / 2 / sqrt(fabs(sin(t)))) / 3;
		else return 2 * (-Q * sin(t) - sqrt(Q) * cos(t) / sqrt(fabs(sin(t)))) / 3;
	}
} y;

class S : public Function<double>{
public:
	virtual double operator () (const double& t) const {
		return sqrt(x.d(t) * x.d(t) + y.d(t) * y.d(t));
	}
} s;

// (l, r) 上的辛普森积分
double IS(const Function<double>& f, const double& l, const double& r) {
	return (r-l) * (f(l) + 4 * f((l+r)/2) + f(r)) / 6;
}
// 自适应辛普森积分
double Simpson(const Function<double>& f, const double& l, const double& r) {
	double mid = (l+r) / 2, s = IS(f, l, r), t = IS(f, l, mid) + IS(f, mid, r);
	if (fabs(s-t) < eps) return s;
	return Simpson(f, l, mid) + Simpson(f, mid, r);
}

// 求出曲线的 n 等分点
vector<double> Get_Knots(int n) {
	double I = Simpson(s, dt, PI - dt) * 2;
	vector<double> x(n+1), y(n+1);
	x[0] = dt;
	// 在曲线上按长度均匀取插值结点
	for (int i = 1; i <= n/2-1; ++ i) {
		double l = x[i-1], r = 2*PI;
		while(r-l>=1e-8){
			double mid = (l+r) / 2;
			double Im = Simpson(s, x[i-1], mid);
			if (Im < I/n) l = mid;
			else r = mid;
		}
		x[i] = l;
	}
	x[n/2] = PI - dt;
	for (int i = 0; i <= n/2-1; ++ i) {
		x[n-i] = 2*PI - x[i];
	}
	return x;
}

// 用线性样条拟合曲线
// file 曲线数据的文件名
void Linear_heart_plot(int n, int m, const char* file) {
	vector<double> t = Get_Knots(n), f(n+1), g(n+1);
	for (int i = 0; i <= n; ++ i) f[i] = x(t[i]);
	for (int i = 0; i <= n; ++ i) g[i] = y(t[i]);
	LinearSpline<double> X = LinearSpline<double>(t, f);
	LinearSpline<double> Y = LinearSpline<double>(t, g);
	ofstream out(file);
	out << "x,y" << '\n';
	for (int i = 0; i <= m; ++ i) out << X.GetValue(dt+(2*PI-2*dt)*i/m) << ',' << Y.GetValue(dt+(2*PI-2*dt)*i/m) << '\n';
}

// 用三次样条拟合曲线（这里只实现了自然三次样条）
void Cubic_heart_plot(int n, int m, const string& mode, const char* file) {
	vector<double> t = Get_Knots(n), f(n+1), g(n+1);
	for (int i = 0; i <= n; ++ i) f[i] = x(t[i]);
	for (int i = 0; i <= n; ++ i) g[i] = y(t[i]);
	double m0 = 0, mn = 0;
	CubicSpline<double> X = CubicSpline<double>(t, f, mode, m0, mn);
	CubicSpline<double> Y = CubicSpline<double>(t, g, mode, m0, mn);
	ofstream out(file);
	out << "x,y" << '\n';
	for (int i = 0; i <= m; ++ i) out << X.GetValue(dt+(2*PI-2*dt)*i/m) << ',' << Y.GetValue(dt+(2*PI-2*dt)*i/m) << '\n';
}

int main() {
	Linear_heart_plot(10, 50000, "e1.csv");
	Linear_heart_plot(40, 50000, "e2.csv");
	Linear_heart_plot(160, 50000, "e3.csv");
	Cubic_heart_plot(10, 50000, "Natural", "e4.csv");
	Cubic_heart_plot(40, 50000, "Natural", "e5.csv");
	Cubic_heart_plot(160, 50000, "Natural", "e6.csv");
}