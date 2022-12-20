#include <bits/stdc++.h>
#include "Spline.h"
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

void Linear_heart_plot(int n, int m, const char* file) {
	LinearSpline<double> X = LinearSplineInterpolation<double>(x, 0, 2*PI, n);
	LinearSpline<double> Y = LinearSplineInterpolation<double>(y, 0, 2*PI, n);
	ofstream out(file);
	out << "x,y" << '\n';
	for (int i = 0; i <= m; ++ i) out << X.GetValue(2*PI*i/m) << ',' << Y.GetValue(2*PI*i/m) << '\n';
}

void Cubic_heart_plot(int n, int m, const string& mode, const char* file) {
	ofstream out(file);
	out << "x,y" << '\n';
	if (mode == "Natural") {
		CubicSpline<double> X = CubicSplineInterpolation<double>(x, 0, 2*PI, n, mode);
		CubicSpline<double> Y = CubicSplineInterpolation<double>(y, 0, 2*PI, n, mode);
		for (int i = 0; i <= m; ++ i) out << X.GetValue(2*PI*i/m) << ',' << Y.GetValue(2*PI*i/m) << '\n';
	}
	else {
		CubicSpline<double> X = CubicSplineInterpolation<double>(x, PI*0.2, PI*2.2, n, mode);
		CubicSpline<double> Y = CubicSplineInterpolation<double>(y, PI*0.2, PI*2.2, n, mode);
		for (int i = 0; i <= m; ++ i) out << X.GetValue(PI*0.2+2*PI*i/m) << ',' << Y.GetValue(PI*0.2+2*PI*i/m) << '\n';
	}
}

int main() {
	Linear_heart_plot(10, 5000, "e1.csv");
	Linear_heart_plot(40, 5000, "e2.csv");
	Linear_heart_plot(160, 5000, "e3.csv");
	Cubic_heart_plot(10, 5000, "Natural", "e4.csv");
	Cubic_heart_plot(40, 5000, "Natural", "e5.csv");
	Cubic_heart_plot(160, 5000, "Natural", "e6.csv");
	Cubic_heart_plot(10, 5000, "Complete", "e7.csv");
	Cubic_heart_plot(40, 5000, "Complete", "e8.csv");
	Cubic_heart_plot(160, 5000, "Complete", "e9.csv");
}