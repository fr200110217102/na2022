#ifndef GEO2D
#define GEO2D
#include <bits/stdc++.h>
#include "function.h"
using namespace std;

inline int sgn(double x) {
	if (fabs(x) < 1e-12) return 0;
	return x>0 ? 1 : -1;
}

struct vec {
	double x,y;
	vec() {}
	vec(double x,double y):x(x),y(y){}
	vec operator + (const vec& p) const {
		return vec(x+p.x, y+p.y);
	}
	vec operator - (const vec& p) const {
		return vec(x-p.x, y-p.y);
	}
	double norm() const {
		return sqrt(x*x + y*y);
	}
	bool operator < (const vec& p) const {
		return sgn(x-p.x) == -1 || (sgn(x-p.x) == 0 && sgn(y-p.y) == -1);
	}
	bool operator == (const vec& p) const {
		return sgn(x-p.x) == 0 && sgn(y-p.y) == 0;
	}
};
typedef vec pnt;

inline double dis(const pnt& p, const pnt& q) {
	return (p-q).norm();
}

struct circle {
	pnt o;
	double r;
	circle(double x=0, double y=0, double r=0) : o(x,y), r(r) {}
	double Y1(const double& x) const {return o.y - sqrt(max(0.0, r*r - (x-o.x)*(x-o.x)));}
	double Y2(const double& x) const {return o.y + sqrt(max(0.0, r*r - (x-o.x)*(x-o.x)));}
	double X1(const double& y) const {return o.x - sqrt(max(0.0, r*r - (y-o.y)*(y-o.y)));}
	double X2(const double& y) const {return o.x + sqrt(max(0.0, r*r - (y-o.y)*(y-o.y)));}
	double arg(const pnt& p) const {return atan2(p.x-o.x, p.y-o.y);}
};

inline int pnt_to_circle(const pnt &p, const circle &c) {
	return sgn(dis(p, c.o) - c.r);
}
#endif