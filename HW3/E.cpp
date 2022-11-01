#include <bits/stdc++.h>
#include "BSpline.h"
#include "../HW1/function.h"
using namespace std;

const double PI = acos(-1), Q = sqrt(3);

class X : public Function<double> {
public:
	virtual double operator () (const double& t) const {
		return sin(t);
	}
} x;

class Y : public Function<double>{
public:
	virtual double operator () (const double& t) const {
		return 2 * (Q * cos(t) + sqrt(Q * fabs(sin(t)))) / 3;
	}
} y;

void heart_plot(int n, int m, const char* file) {
	vector <double> f;
	f.resize(n/2);
	for (int i = 0; i <= n/2 - 1; ++ i) f[i] = x((i+0.5)*2*PI/n);
	QuadraticBSpline<double> X1 = QuadraticBSpline<double>(0, n/2, f, 0, 0);
	for (int i = 0; i <= n/2 - 1; ++ i) f[i] = y((i+0.5)*2*PI/n);
	QuadraticBSpline<double> Y1 = QuadraticBSpline<double>(0, n/2, f, 2/Q, -2/Q);
	for (int i = 0; i <= n/2 - 1; ++ i) f[i] = x((i+n/2+0.5)*2*PI/n);
	QuadraticBSpline<double> X2 = QuadraticBSpline<double>(0, n/2, f, 0, 0);
	for (int i = 0; i <= n/2 - 1; ++ i) f[i] = y((i+n/2+0.5)*2*PI/n);
	QuadraticBSpline<double> Y2 = QuadraticBSpline<double>(0, n/2, f, -2/Q, 2/Q);
	ofstream out1(file);
	out1 << "x,y" << '\n';
	for (int i = 0; i < m/2; ++ i) out1 << X1.GetValue(1.0*i*n/m) << ',' << Y1.GetValue(1.0*i*n/m) << '\n';
	for (int i = 0; i < m/2; ++ i) out1 << X2.GetValue(1.0*i*n/m) << ',' << Y2.GetValue(1.0*i*n/m) << '\n';

}

int main() {
	heart_plot(10, 5000, "e1.csv");
	heart_plot(40, 5000, "e2.csv");
	heart_plot(160, 5000, "e3.csv");
}