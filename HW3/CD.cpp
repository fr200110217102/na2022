#include <bits/stdc++.h>
#include "BSpline.h"
#include "../HW1/function.h"
using namespace std;

class F : public Function<double> {
public:
	virtual double operator () (const double& x) const {
		return 1 / (1 + x*x);
	}
	virtual double d (const double& x) const {
		return -2*x / (1 + x*x) / (1 + x*x);
	}
} f;

int main() {
	CubicBSpline<double> C1 = CubicBSplineInterpolation<double>(f, -5, 10);
	ofstream out1("c1.csv");
	out1 << "x,y" << '\n';
	for(int i = -4999; i <= 4999; ++ i) out1 << i * 0.001 << ',' << C1.GetValue(i * 0.001) << '\n';
	
	QuadraticBSpline<double> C2 = QuadraticBSplineInterpolation<double>(f, -5, 10);
	ofstream out2("c2.csv");
	out2 << "x,y" << '\n';
	for(int i = -4999; i <= 4999; ++ i) out2 << i * 0.001 << ',' << C2.GetValue(i * 0.001) << '\n';

	ofstream outd("d.csv");
	vector<double> a = {-3.5,-3,-0.5,0,0.5,3,3.5};
	outd << "a_i,ES3,ES2" << '\n';
	for(int i = 0; i < 7; ++ i) outd << a[i] << ',' << fabs(C1.GetValue(a[i]) - f(a[i])) << ',' << fabs(C2.GetValue(a[i]) - f(a[i])) << '\n';
}