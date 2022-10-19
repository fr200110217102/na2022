#include <bits/stdc++.h>
#include "function.h"
#include "EquationSolver.h"
using namespace std;

const double PI = acos(-1);

class F1 : public Function <double> {
    virtual double operator () (const double& x) const {
        return 1.0 / x - tan(x);
    }
}f1;

class F2 : public Function <double> {
    virtual double operator () (const double& x) const {
        return 1.0 / x - pow(2, x);
    }
}f2;

class F3 : public Function <double> {
    virtual double operator () (const double& x) const {
        return pow(2, -x) + exp(x) + 2 * cos(x) - 6;
    }
}f3;

class F4 : public Function <double> {
    virtual double operator () (const double& x) const {
        return (((x + 4) * x + 3) * x + 5) / (((2 * x - 9) * x + 18) * x - 2);
    }
}f4;

class G : public Function <double> {
    virtual double operator () (const double& x) const {
        return x - tan(x);
    }
    virtual double d(const double& x) const {
        double t = cos(x);
        return 1 - 1.0 / (t * t);
    }
}g;

class H1 : public Function <double> {
    virtual double operator () (const double& x) const {
        return sin(x / 2) - 1;
    }
}h1;

class H2 : public Function <double> {
    virtual double operator () (const double& x) const {
        return exp(x) - tan(x);
    }
}h2;

class H3 : public Function <double> {
    virtual double operator () (const double& x) const {
        return ((x - 12) * x + 3) * x + 1;
    }
}h3;

class P : public Function <double> {
    virtual double operator () (const double& h) const {
        return L * (PI/2 * r * r - r * r * asin(h / r) - h * sqrt(r * r - h * h)) - V;
    }
    virtual double d(const double& h) const {
        return L * (-2 * sqrt(r * r - h * h));
    }
private:
    double L, r, V;
public:
    P(double L, double r, double V) : L(L), r(r), V(V) {}
};

class Q : public Function <double> {
    virtual double operator () (const double& a) const {
        double _a = a * PI / 180, s = sin(_a), c = cos(_a);
        return A * s * c + B * s * s - C * c - E * s;
    }
    virtual double d(const double& a) const {
        double _a = a * PI / 180;
        return (A * cos(2 * _a) + B * sin(2 * _a) + C * sin(_a) - E * cos(_a)) * (PI/180);
    }
private:
    double A, B, C, E;
public:
    Q(const double& l, const double& h, const double& D, const double& b1) {
        double _b = b1 * PI / 180;
        A = l * sin(_b);
        B = l * cos(_b);
        C = (h + 0.5 * D) * sin(_b) - 0.5 * D * tan(_b);
        E = (h + 0.5 * D) * cos(_b) - 0.5 * D;
    }
};

int main() {
    cout << "2(1)\n" << Bisection<double>(f1, 0.0, PI/2).solve() << endl;
    cout << "2(2)\n" << Bisection<double>(f2, 0.0, 1.0).solve() << endl;
    cout << "2(3)\n" << Bisection<double>(f3, 1.0, 3.0).solve() << endl;
    cout << "2(4)\n" << Bisection<double>(f4, 0.0, 4.0).solve() << endl;

    cout << "3(1)\n" << Newton<double>(g, 4.5).solve() << endl;
    cout << "3(2)\n" << Newton<double>(g, 7.7).solve() << endl;

    cout << "4(1)\n" << Secant<double>(h1, 0.0, PI/2).solve() << endl;
    cout << "4(2)\n" << Secant<double>(h2, 1.0, 1.4).solve() << endl;
    cout << "4(3)\n" << Secant<double>(h3, 0.0, -0.5).solve() << endl;

    cout << "5\n";
    P p(10, 1, 12.4);
    cout << Bisection<double>(p, 0, 1, 20, 0.001).solve() << endl;
    cout << Newton<double>(p, 0.5).solve() << endl;
    cout << Secant<double>(p, 0, 1, 20, 0.001).solve() << endl;
    
    cout << "6\n";
    Q q(89, 49, 55, 11.5);
    cout << Newton<double>(q, 33).solve() << endl;
    q = Q(89, 49, 30, 11.5);
    cout << Newton<double>(q, 33).solve() << endl;

    cout << Secant<double>(q, 30, 45).solve() << endl;
    cout << Secant<double>(q, 60, 90).solve() << endl;
    cout << Secant<double>(q, 90, 180).solve() << endl;
    cout << Secant<double>(q, 180, 360).solve() << endl;
}