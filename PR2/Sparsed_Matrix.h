#ifndef SPARSEDMATRIX
#define SPARSEDMATRIX
#include <bits/stdc++.h>
#include "Matrix.h"
using namespace std;

template<class type>
struct Sparsed_Matrix {
	int n, m;
	vector<map<int,type>> ele;
	Sparsed_Matrix(int n, int m): n(n), m(m) {ele.resize(n);}
	map<int,type>& operator[](int i) {
		return ele[i];
	}
	const map<int,type>& operator[](int i) const {
		return ele[i];
	}
	Colvec<type> operator*(const Colvec<type>& v)const{
		Colvec<type> res(n);
		for (int i = 0; i < n; ++ i)
			for (auto& [j, x] : ele[i])
				res[i] += v[j] * x;
		return res;
	}
	operator Matrix<type>() const {
		Matrix<type> res(n, m);
		for (int i = 0; i < n; ++ i)
			for (auto& [j, x] : ele[i])
				res[i][j] = x;
		return res;
	}
};
#endif