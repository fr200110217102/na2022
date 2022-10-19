#include <bits/stdc++.h>
#include "CubicSpine.h"
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
	CubicSpine<double> S1 = CubicSpineInterpolation<double>(f, -1.0, 1.0, 5);
	CubicSpine<double> S2 = CubicSpineInterpolation<double>(f, -1.0, 1.0, 10);
	CubicSpine<double> S3 = CubicSpineInterpolation<double>(f, -1.0, 1.0, 20);
	CubicSpine<double> S4 = CubicSpineInterpolation<double>(f, -1.0, 1.0, 40);
	CubicSpine<double> S5 = CubicSpineInterpolation<double>(f, -1.0, 1.0, 80);
	
	ofstream out1("s1.csv");
	out1 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out1 << i * 0.001 << ',' << S1.GetValue(i * 0.001) << '\n';
	ofstream out2("s2.csv");
	out2 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out2 << i * 0.001 << ',' << S2.GetValue(i * 0.001) << '\n';
	ofstream out3("s3.csv");
	out3 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out3 << i * 0.001 << ',' << S3.GetValue(i * 0.001) << '\n';
	ofstream out4("s4.csv");
	out4 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out4 << i * 0.001 << ',' << S4.GetValue(i * 0.001) << '\n';
	ofstream out5("s5.csv");
	out5 << "x,y" << '\n';
	for(int i = -999; i <= 999; ++ i) out5 << i * 0.001 << ',' << S5.GetValue(i * 0.001) << '\n';
}