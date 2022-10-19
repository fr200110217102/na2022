#include <bits/stdc++.h>
#include "../HW1/function.h"
#include "PolynomialInterpolation.h"
using namespace std;

class F1 : public Function <double> {
public:
    virtual double operator () (const double &x) const {
        return 1 / (x * x + 1);
    }
} f1;

class F2 : public Function <double> {
public:
    virtual double operator () (const double &x) const {
        return 1 / (25 * x * x + 1);
    }
} f2;

int main () {
    Newton_Interpolation<double> ip1_2 = Grid_Interpolation(f1, -5, 5, 2);
    Newton_Interpolation<double> ip1_4 = Grid_Interpolation(f1, -5, 5, 4);
    Newton_Interpolation<double> ip1_6 = Grid_Interpolation(f1, -5, 5, 6);
    Newton_Interpolation<double> ip1_8 = Grid_Interpolation(f1, -5, 5, 8);
    Polynomial<double> p1_2 = ip1_2.GetPolynomial();
    Polynomial<double> p1_4 = ip1_4.GetPolynomial();
    Polynomial<double> p1_6 = ip1_6.GetPolynomial();
    Polynomial<double> p1_8 = ip1_8.GetPolynomial();
    cout << p1_2 << endl;
    cout << p1_4 << endl;
    cout << p1_6 << endl;
    cout << p1_8 << endl;

    Newton_Interpolation<double> ip2_2 = Chebyshev_Interpolation(f2, -1, 1, 5);
    Newton_Interpolation<double> ip2_4 = Chebyshev_Interpolation(f2, -1, 1, 10);
    Newton_Interpolation<double> ip2_6 = Chebyshev_Interpolation(f2, -1, 1, 15);
    Newton_Interpolation<double> ip2_8 = Chebyshev_Interpolation(f2, -1, 1, 20);
    Polynomial<double> p2_2 = ip2_2.GetPolynomial();
    Polynomial<double> p2_4 = ip2_4.GetPolynomial();
    Polynomial<double> p2_6 = ip2_6.GetPolynomial();
    Polynomial<double> p2_8 = ip2_8.GetPolynomial();
    cout << p2_2 << endl;
    cout << p2_4 << endl;
    cout << p2_6 << endl;
    cout << p2_8 << endl;

    vector <double> x(10), y(10);
    x[0] = x[1] = 0, x[2] = x[3] = 3, x[4] = x[5] = 5, x[6] = x[7] = 8, x[8] = x[9] = 13;
    y[0] = 0, y[2] = 225, y[4] = 383, y[6] = 623, y[8] = 993;
    y[1] = 75, y[3] = 77, y[5] = 80, y[7] = 74, y[9] = 72;
    Hermite_Interpolation<double> ip3(x, y);
    Polynomial<double> p3 = ip3.GetPolynomial();
    cout << p3 << endl;
    cout << p3.d() << endl;
    cout << p3(10) << ' ' << p3.d(10) << endl;

    x.resize(7);
    vector <double> y1(7), y2(7);
    x[0] = 0, x[1] = 6, x[2] = 10, x[3] = 13, x[4] = 17, x[5] = 20, x[6] = 28;
    y1[0] = 6.67, y1[1] = 17.3, y1[2] = 42.7, y1[3] = 37.3, y1[4] = 30.1, y1[5] = 29.3, y1[6] = 28.7;
    y2[0] = 6.67, y2[1] = 16.1, y2[2] = 18.9, y2[3] = 15.0, y2[4] = 10.6, y2[5] = 9.44, y2[6] = 8.89;
    Newton_Interpolation<double> ip4_1(x, y1);
    Polynomial<double> p4_1 = ip4_1.GetPolynomial();
    Newton_Interpolation<double> ip4_2(x, y2);
    Polynomial<double> p4_2 = ip4_2.GetPolynomial();
    cout << p4_1 << endl;
    cout << p4_2 << endl;
    cout << p4_1(43) << ' ' << p4_2(43) << endl;
}