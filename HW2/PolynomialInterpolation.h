#include <bits/stdc++.h>
#include "../HW1/function.h"
using namespace std;
const double PI = 3.14159265358979323846;

template <class type> 
class Polynomial : public vector <type>, public Function <type> {
public :
    Polynomial (const int & n = -1, const type & x0 = 0, const type & x1 = 0) {
        this -> resize(n+1);
        if (n >= 0) this -> at(0) = x0;
        if (n >= 1) this -> at(1) = x1;
    }
    int time() const {
        return this -> size() - 1;
    }
    Polynomial <type> operator + (const Polynomial <type>& b) const {
        int n = max(time(), b.time()), m = min(time(), b.time());
        Polynomial <type> c(n);
        for (int i = 0; i <= m; ++ i) c[i] = this -> at(i) + b[i];
        if (n == time())
            for (int i = m+1; i <= n; ++ i) c[i] = this -> at(i);
        else
            for (int i = m+1; i <= n; ++ i) c[i] = b[i];
        return c;
    }
    Polynomial <type> operator += (const Polynomial <type>& b) {
        return *this = *this + b;
    }
    Polynomial <type> operator * (const type& b) const {
        int n = time();
        Polynomial <type> c(n);
        for (int i = 0; i <= n; ++ i) c[i] = this -> at(i) * b;
        return c;
    }
    Polynomial <type> operator *= (const type & b) {
        return *this = *this * b;
    }
    Polynomial <type> operator * (const Polynomial <type>& b) const {
        int n = time(), m = b.time();
        Polynomial <type> c(n+m);
        for (int i = 0; i <= n; ++ i)
            for (int j = 0; j <= m; ++ j)
                c[i+j] += this -> at(i) * b[j];
        return c;
    }
    Polynomial <type> operator *= (const Polynomial <type>& b) {
        return *this = *this * b;
    }
    virtual type operator ()(const type & x) const {
        int n = time();
        type res = 0;
        for (int i = n; i >= 0; -- i)
            res = res * x + this -> at(i);
        return res;
    }
    Polynomial <type> d() {
        int n = time();
        Polynomial <type> c(n-1);
        for (int i = n; i > 0; -- i)
            c[i-1] = this -> at(i) * i;
        return c;
    }
    virtual type d(const type & x) const {
        int n = time();
        type res = 0;
        for (int i = n; i > 0; -- i)
            res = res * x + this -> at(i) * i;
        return res;
    }
};

template <class type>
Polynomial <type> operator * (const type & x, const Polynomial <type> & p) {
    return p * x;
}

template <class type>
ostream & operator << (ostream &out, const Polynomial <type> &p) {
    int n = p.time();
    if (n == -1) {cout << 0; return out;}
    for (int i = 0; i <= n; ++ i) {
        out << p[i];
        if (i == 1) out << "*x";
        else if (i >= 2) out << "*x**" << i;
        if (i < n) {
            if (p[i+1] >= 0) out << " +";
            else cout << " ";
        }
    }
    return out;
}

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