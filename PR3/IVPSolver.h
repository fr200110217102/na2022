#ifndef __IVPSOLVER__
#define __IVPSOLVER__
#include <bits/stdc++.h>
#include "Function.h"
using namespace std;

class IVPSolver {
protected:
	string method_name;
	vector<Colvec<double>> u;
public:
	virtual vector<Colvec<double>> Solve(const Function_nd<double>& f, const Colvec<double>& u0, double T, int n, int s = 0, bool omit = 0) = 0;
};
#endif