#ifndef BVPSOLVER
#define BVPSOLVER
#include <bits/stdc++.h>
#include "function.h"
using namespace std;

class BVPSolver {
protected:
	virtual void Matrix_Generator() = 0;
public:
	virtual void Solve() = 0;
	virtual double getval(int, int) = 0;
};
#endif