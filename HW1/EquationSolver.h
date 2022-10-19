#ifndef EQUATIONSOLVER
#define EQUATIONSOLVER
#include <iostream>
#include <algorithm>
#include <limits>
#include "function.h"

const double eps = 1e-12;

template <class type>
class EquationSolver{
protected:
    virtual type solve() = 0;
};

template <class type>
class Bisection : public EquationSolver <type> {
private:
    const Function<type> &f;
    type a, b, delta;
    int M;
public:
    Bisection(const Function<type> &f, const type &a, const type &b, const int &M = 100, const type &delta = 1e-6) :
        f(f), a(a), b(b), delta(delta), M(M) {}
    virtual type solve() {
        if (f(a) * f(b) > eps) throw "Invalid Interval!";
        type h = b - a, u = f(a), c, w, x = a;
        int k = 1;
        while (k <= M) {
            h /= 2, c = x + h, w = f(c);
            if (fabs(h) < delta || fabs(w) < eps) break;
            else if (w * u > 0) x = c;
            ++ k;
        }
        if (k > M) std::cout << "Time Limit Exceeded!" << std::endl;
        std::cerr << "Bisection : times = " << k << ", " << "delta = " << h << std::endl;
        return c;
    }
};

template <class type>
class Newton : public EquationSolver <type> {
private:
    const Function<type> &f;
    type x0;
    int M;
public:
    Newton(const Function<type> &f, const type &x0, const int& M = 10) :
        f(f), x0(x0), M(M){}
    virtual type solve() {
        type x = x0, u;
        int k = 1;
        while (k <= M) {
            u = f(x);
            if (fabs(u) < eps) break;
            x -= u / f.d(x);
            ++ k;
        }
        if (k > M) std::cout << "Time Limit Exceeded!" << std::endl;
        std::cerr << "Newton : times = " << k << std::endl;
        return x;
    }
};

template <class type>
class Secant : public EquationSolver <type> {
private:
    const Function<type> &f;
    type a, b, delta;
    int M;
public:
    Secant<type>(const Function<type> &f, const type &a, const type &b, const int& M = 30, const type& delta = 1e-6) :
        f(f), a(a), b(b), delta(delta), M(M) {}
    virtual type solve() {
        type x0 = a, x1 = b, u = f(x1), v = f(x0), s;
        int k = 2;
        while (k <= M) {
            if (fabs(u) > fabs(v)) std::swap(x0, x1), std::swap(u, v);
            s = (x1 - x0) / (u - v);
            x0 = x1, v = u;
            x1 -= u * s, u = f(x1);
            if (fabs(x0 - x1) < delta || fabs(u) < eps) break;
            ++ k;
        }
        if (k > M) std::cout << "Time Limit Exceeded!" << std::endl;
        std::cerr << "Secant : times = " << k << ", " << "delta = " << fabs(x0 - x1) << std::endl;
        return x1;
    }
};
#endif