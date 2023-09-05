#include <bits/stdc++.h>
#include "Function.h"
#include "Matrix.h"
#include "MOL.h"
using namespace std;

class LeapFrog : public MOL_for_adv {
protected:
	virtual void step(const Function<double>& f, const int& n) {
		double P, Q, L, R;
		if (n == 0) {
			for (int i = 0; i <= m; ++ i) {
				P = u[n][i], L = u[n][i==0?m:i-1], R = u[n][i==m?0:i+1];
				u[n+1][i] = P - mu*.5 * (R - L);
			}
		}
		else {
			for (int i = 0; i <= m; ++ i) {
				Q = u[n-1][i], L = u[n][i==0?m:i-1], R = u[n][i==m?0:i+1];
				u[n+1][i] = Q - mu * (R - L);
			}
		}
	}
};

class Lax_Friedrichs : public MOL_for_adv {
protected:
	virtual void step(const Function<double>& f, const int& n) {
		double P, L, R;
		for (int i = 0; i <= m; ++ i) {
			P = u[n][i], L = u[n][i==0?m:i-1], R = u[n][i==m?0:i+1];
			u[n+1][i] = .5 * (L + R) - mu*.5 * (R - L);
		}
	}
};

class Lax_Wendroff : public MOL_for_adv {
protected:
	virtual void step(const Function<double>& f, const int& n) {
		double P, L, R;
		for (int i = 0; i <= m; ++ i) {
			P = u[n][i], L = u[n][i==0?m:i-1], R = u[n][i==m?0:i+1];
			u[n+1][i] = P - mu*.5 * (R - L) + mu*mu*.5 * (L - 2*P + R);
		}
	}
};

class Upwind : public MOL_for_adv {
protected:
	virtual void step(const Function<double>& f, const int& n) {
		double P, L, R;
		if (mu >= 0) {
			for (int i = 0; i <= m; ++ i) {
				P = u[n][i], L = u[n][i==0?m:i-1];
				u[n+1][i] = P - mu * (P - L);
			}
		}
		else {
			for (int i = 0; i <= m; ++ i) {
				P = u[n][i], R = u[n][i==m?0:i+1];
				u[n+1][i] = P - mu * (R - P);
			}
		}
	}
};

class Beam_Warming : public MOL_for_adv {
protected:
	virtual void step(const Function<double>& f, const int& n) {
		double P, L, LL, R, RR;
		if (mu >= 0) {
			for (int i = 0; i <= m; ++ i) {
				P = u[n][i], L = u[n][i==0?m:i-1], LL = u[n][i<=1?m+i-2:i-2];
				u[n+1][i] = P - mu*.5 * (3*P - 4*L + LL) + mu*mu*.5 * (P - 2*L + LL);
			}
		}
		else {
			for (int i = 0; i <= m; ++ i) {
				P = u[n][i], R = u[n][i==m?0:i+1], RR = u[n][i>=m-1?i+2-m:i+2];
				u[n+1][i] = P - mu*.5 * (3*P - 4*R + RR) + mu*mu*.5 * (P - 2*R + RR);
			}
		}
	}
};