#include <bits/stdc++.h>
#include "CubicSpline.h"
#include "../HW1/function.h"
using namespace std;

class F : public Function<double> {
	virtual double operator () (const double& x) const {
		return 1 / (1 + 25*x*x);
	}
	virtual double d (const double& x) const {
		return -50 * x / (1 + 25*x*x) / (1 + 25*x*x);
	}
} f;

int main() {
	CubicSpline<double> A1 = CubicSplineInterpolation<double>(f, -1.0, 1.0, 5);
	CubicSpline<double> A2 = CubicSplineInterpolation<double>(f, -1.0, 1.0, 10);
	CubicSpline<double> A3 = CubicSplineInterpolation<double>(f, -1.0, 1.0, 20);
	CubicSpline<double> A4 = CubicSplineInterpolation<double>(f, -1.0, 1.0, 40);
	CubicSpline<double> A5 = CubicSplineInterpolation<double>(f, -1.0, 1.0, 80);
	
	ofstream out1("a1.csv");
	out1 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out1 << i * 0.001 << ',' << A1.GetValue(i * 0.001) << '\n';
	ofstream out2("a2.csv");
	out2 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out2 << i * 0.001 << ',' << A2.GetValue(i * 0.001) << '\n';
	ofstream out3("a3.csv");
	out3 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out3 << i * 0.001 << ',' << A3.GetValue(i * 0.001) << '\n';
	ofstream out4("a4.csv");
	out4 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out4 << i * 0.001 << ',' << A4.GetValue(i * 0.001) << '\n';
	ofstream out5("a5.csv");
	out5 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out5 << i * 0.001 << ',' << A5.GetValue(i * 0.001) << '\n';
}