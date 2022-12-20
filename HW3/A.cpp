#include <bits/stdc++.h>
#include "Spline.h"
#include "../HW1/function.h"
using namespace std;

class F : public Function<double> {
	virtual double operator () (const double& x) const {
		return 1 / (1 + 25*x*x);
	}
	virtual double d (const double& x, const int& k = 1) const {
		if (k == 1) return -50 * x / (1 + 25*x*x) / (1 + 25*x*x);
		else if (k == 2) return (-50 + 3750*x*x) / (1 + 25*x*x) / (1 + 25*x*x) / (1 + 25*x*x);
		else throw 0;
	}
} f;

int main() {
	LinearSpline<double> A1 = LinearSplineInterpolation<double>(f, -1.0, 1.0, 5);
	LinearSpline<double> A2 = LinearSplineInterpolation<double>(f, -1.0, 1.0, 10);
	LinearSpline<double> A3 = LinearSplineInterpolation<double>(f, -1.0, 1.0, 20);
	LinearSpline<double> A4 = LinearSplineInterpolation<double>(f, -1.0, 1.0, 40);
	LinearSpline<double> A5 = LinearSplineInterpolation<double>(f, -1.0, 1.0, 80);

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

	CubicSpline<double> A6 = CubicSplineInterpolation<double>(f, -1.0, 1.0, 5, "Complete");
	CubicSpline<double> A7 = CubicSplineInterpolation<double>(f, -1.0, 1.0, 10, "Complete");
	CubicSpline<double> A8 = CubicSplineInterpolation<double>(f, -1.0, 1.0, 20, "Complete");
	CubicSpline<double> A9 = CubicSplineInterpolation<double>(f, -1.0, 1.0, 40, "Complete");
	CubicSpline<double> A10 = CubicSplineInterpolation<double>(f, -1.0, 1.0, 80, "Complete");
	
	ofstream out6("a6.csv");
	out6 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out6 << i * 0.001 << ',' << A6.GetValue(i * 0.001) << '\n';
	ofstream out7("a7.csv");
	out7 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out7 << i * 0.001 << ',' << A7.GetValue(i * 0.001) << '\n';
	ofstream out8("a8.csv");
	out8 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out8 << i * 0.001 << ',' << A8.GetValue(i * 0.001) << '\n';
	ofstream out9("a9.csv");
	out9 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out9 << i * 0.001 << ',' << A9.GetValue(i * 0.001) << '\n';
	ofstream out10("a10.csv");
	out10 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out10 << i * 0.001 << ',' << A10.GetValue(i * 0.001) << '\n';

	CubicSpline<double> A11 = CubicSplineInterpolation<double>(f, -1.0, 1.0, 5, "Natural");
	CubicSpline<double> A12 = CubicSplineInterpolation<double>(f, -1.0, 1.0, 10, "Natural");
	CubicSpline<double> A13 = CubicSplineInterpolation<double>(f, -1.0, 1.0, 20, "Natural");
	CubicSpline<double> A14 = CubicSplineInterpolation<double>(f, -1.0, 1.0, 40, "Natural");
	CubicSpline<double> A15 = CubicSplineInterpolation<double>(f, -1.0, 1.0, 80, "Natural");
	
	ofstream out11("a11.csv");
	out11 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out11 << i * 0.001 << ',' << A11.GetValue(i * 0.001) << '\n';
	ofstream out12("a12.csv");
	out12 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out12 << i * 0.001 << ',' << A12.GetValue(i * 0.001) << '\n';
	ofstream out13("a13.csv");
	out13 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out13 << i * 0.001 << ',' << A13.GetValue(i * 0.001) << '\n';
	ofstream out14("a14.csv");
	out14 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out14 << i * 0.001 << ',' << A14.GetValue(i * 0.001) << '\n';
	ofstream out15("a15.csv");
	out15 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out15 << i * 0.001 << ',' << A15.GetValue(i * 0.001) << '\n';

	CubicSpline<double> A16 = CubicSplineInterpolation<double>(f, -1.0, 1.0, 5, "Specified_Second_Derivatives");
	CubicSpline<double> A17 = CubicSplineInterpolation<double>(f, -1.0, 1.0, 10, "Specified_Second_Derivatives");
	CubicSpline<double> A18 = CubicSplineInterpolation<double>(f, -1.0, 1.0, 20, "Specified_Second_Derivatives");
	CubicSpline<double> A19 = CubicSplineInterpolation<double>(f, -1.0, 1.0, 40, "Specified_Second_Derivatives");
	CubicSpline<double> A20 = CubicSplineInterpolation<double>(f, -1.0, 1.0, 80, "Specified_Second_Derivatives");

	ofstream out16("a16.csv");
	out16 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out16 << i * 0.001 << ',' << A16.GetValue(i * 0.001) << '\n';
	ofstream out17("a17.csv");
	out17 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out17 << i * 0.001 << ',' << A17.GetValue(i * 0.001) << '\n';
	ofstream out18("a18.csv");
	out18 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out18 << i * 0.001 << ',' << A18.GetValue(i * 0.001) << '\n';
	ofstream out19("a19.csv");
	out19 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out19 << i * 0.001 << ',' << A19.GetValue(i * 0.001) << '\n';
	ofstream out20("a20.csv");
	out20 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out20 << i * 0.001 << ',' << A20.GetValue(i * 0.001) << '\n';

}