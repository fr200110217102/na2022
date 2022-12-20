#include <bits/stdc++.h>
#include "BSpline.h"
#include "../HW1/function.h"
using namespace std;

class F : public Function<double> {
public:
	virtual double operator () (const double& x) const {
		return 1 / (1 + x*x);
	}
	virtual double d (const double& x, const int& k) const {
		if (k == 1) return -2*x / (1 + x*x) / (1 + x*x);
		else if (k == 2) return (-2 + 6*x*x) / (1 + x*x) / (1 + x*x) / (1 + x*x);
		else throw 0;
	}
} f;

int main() {
	LinearBSpline<double> CL = LinearBSplineInterpolation(f, -5, 10);
	ofstream outL("cl.csv");
	outL << "x,y" << endl;
	for(int i = -4999; i <= 4999; ++ i) outL << i * 0.001 << ',' << CL.GetValue(i * 0.001) << endl;

	CubicBSpline<double> C1 = CubicBSplineInterpolation<double>(f, -5, 10, "Natural");
	ofstream out1("c1.csv");
	out1 << "x,y" << endl;
	for(int i = -4999; i <= 4999; ++ i) out1 << i * 0.001 << ',' << C1.GetValue(i * 0.001) << endl;

	CubicBSpline<double> C2 = CubicBSplineInterpolation<double>(f, -5, 10, "Complete");
	ofstream out2("c2.csv");
	out2 << "x,y" << endl;
	for(int i = -4999; i <= 4999; ++ i) out2 << i * 0.001 << ',' << C2.GetValue(i * 0.001) << endl;

	CubicBSpline<double> C3 = CubicBSplineInterpolation<double>(f, -5, 10, "Specified_Second_Derivatives");
	ofstream out3("c3.csv");
	out3 << "x,y" << endl;
	for(int i = -4999; i <= 4999; ++ i) out3 << i * 0.001 << ',' << C3.GetValue(i * 0.001) << endl;
	
	QuadraticBSpline<double> CQ = QuadraticBSplineInterpolation<double>(f, -5, 10);
	ofstream outQ("cq.csv");
	outQ << "x,y" << endl;
	for(int i = -4999; i <= 4999; ++ i) outQ << i * 0.001 << ',' << CQ.GetValue(i * 0.001) << endl;

	ofstream outd("d.csv");
	vector<double> a = {-3.5,-3,-0.5,0,0.5,3,3.5};
	outd << "a_i,ESL,ES1,ES2,ES3,ESQ" << endl;
	for(int i = 0; i < 7; ++ i) outd << a[i] << ',' << fabs(CL.GetValue(a[i]) - f(a[i])) << ',' << fabs(C1.GetValue(a[i]) - f(a[i])) << ',' << fabs(C2.GetValue(a[i]) - f(a[i])) << ',' << fabs(C3.GetValue(a[i]) - f(a[i])) << ',' << fabs(CQ.GetValue(a[i]) - f(a[i])) << endl;
}