#include <bits/stdc++.h>
#include "Matrix.h"
using namespace std;

#ifndef FUNCTION
#define FUNCTION

template <class type>
class Function {
public:
	virtual type operator()(const type& x) const = 0;
	virtual type d(const type& x, const int& i = 1) const = 0;
};

template <class type>
class _0 : public Function<type>{
public:	
	virtual type operator ()(const type& x) const {return 0;}
	virtual type d(const type& x, const int& k = 1) const {return 0;}
};

template <class type>
class Function_nd {
public:
	int n;
	virtual Colvec<type> operator()(const Colvec<type>& x) const = 0;
};

#endif