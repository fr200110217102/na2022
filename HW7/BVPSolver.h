#ifndef BVPSOLVER
#define BVPSOLVER
#include <bits/stdc++.h>
#include "function.h"
using namespace std;

class BVPSolver {
protected:
	virtual void ID_Generator() = 0;
	virtual void Matrix_Generator() = 0;
public:
	virtual void Solve() = 0;
	virtual void Summary(Function_2D<double>&, bool detail) = 0;
};
#endif