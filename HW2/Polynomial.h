#include <bits/stdc++.h>
#include "../HW1/function.h"
using namespace std;

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