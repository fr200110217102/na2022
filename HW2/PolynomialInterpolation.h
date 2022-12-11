#include <bits/stdc++.h>
#include "../HW1/function.h"
#include "Polynomial.h"
using namespace std;
const double PI = 3.14159265358979323846;

template <class type>
class Newton_Interpolation {
private:
    vector <type> x, y;
    vector <vector <type>> f;
public:
    Newton_Interpolation(const vector <type>& x, const vector <type>& y) : x(x), y(y) {
        int n = x.size() - 1;
        f.resize(n+1);
        for (int i = 0; i <= n; ++ i) f[i].resize(n+1);
        for (int i = 0; i <= n; ++ i) f[i][i] = y[i];
        for (int j = 1; j <= n; ++ j) {
            for (int i = 0; i <= n-j; ++ i)
                f[i][i+j] = (f[i+1][i+j] - f[i][i+j-1]) / (x[i+j] - x[i]);
        }
    }
    type GetValue(const type& _x) const {
        int n = x.size() - 1;
        type res = 0, pi = 1;
        for (int i = 0; i <= n; ++ i) {
            res += f[0][i] * pi;
            pi *= _x - x[i];
        }
        return res;
    }
    Polynomial <type> GetPolynomial() const {
        int n = x.size() - 1;
        Polynomial <type> res(0, 0), pi(0, 1);
        for (int i = 0; i <= n; ++ i) {
            res += f[0][i] * pi;
            pi *= Polynomial <type>(1, -x[i], 1);
        }
        return res;
    }
};

template <class type>
class Hermite_Interpolation {
private:
    vector <type> x, y;
    vector <vector <type>> f;
public:
    Hermite_Interpolation(const vector <type>& x, const vector <type>& y) : x(x), y(y) {
        int n = x.size() - 1;
        f.resize(n+1);
        for (int i = 0; i <= n; ++ i) f[i].resize(n+1);
        vector <int> m(n+1);
        for (int i = 0; i <= n; ++ i) {
            if (!i || x[i] != x[i-1]) m[i] = i;
            else m[i] = m[i-1];
        }
        for (int i = 0; i <= n; ++ i) f[i][i] = y[m[i]];
        vector <type> fac(n);
        fac[0] = 1;
        for (int i = 1; i <= n; ++ i) fac[i] = fac[i-1] * i;
        for (int j = 1; j <= n; ++ j)
            for (int i = 0; i <= n-j; ++ i) {
                if (x[i+j] == x[i]) f[i][i+j] = y[m[i]+j] / fac[j];
                else f[i][i+j] = (f[i+1][i+j] - f[i][i+j-1]) / (x[i+j] - x[i]);
            }
    }
    type GetValue(const type &_x) const {
        int n = x.size() - 1;
        type res = 0, pi = 1;
        for (int i = 0; i <= n; ++ i) {
            res += f[0][i] * pi;
            pi *= _x - x[i];
        }
        return res;
    }
    Polynomial <type> GetPolynomial() const {
        int n = x.size() - 1;
        Polynomial <type> res(0, 0), pi(0, 1);
        for (int i = 0; i <= n; ++ i) {
            res += f[0][i] * pi;
            pi *= Polynomial <type>(1, -x[i], 1);
        }
        return res;
    }
};

Newton_Interpolation <double> Grid_Interpolation (const Function <double> & f, const double & l, const double & u, const int & n) {
    vector <double> x(n+1), y(n+1);
    for (int i = 0; i <= n; ++ i) {
        double xi = l + (u - l) * i / n;
        double yi = f(xi);
        x[i] = xi, y[i] = yi;
    }
    return Newton_Interpolation<double>(x, y);
}

Newton_Interpolation <double> Chebyshev_Interpolation (const Function <double> & f, const double & l, const double & u, const int & n) {
    vector <double> x(n+1), y(n+1);
    for (int i = 0; i <= n; ++ i) {
        double xi = (u + l) / 2 + (u - l) / 2 * cos(PI*(2*i+1)/(2*(n+1)));
        double yi = f(xi);
        x[i] = xi, y[i] = yi;
    }
    return Newton_Interpolation<double>(x, y);
}