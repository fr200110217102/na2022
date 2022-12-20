#include <bits/stdc++.h>
using namespace std;

template <class type>
vector<type> Thomas(vector<type>& a, vector<type>& b, vector<type>& c, vector<type> y) {
	int n = a.size();
	if (b.size() != n-1 || c.size() != n-1) throw "Invalid Size!";
	vector<type> p(n), q(n-1);
	p[0] = a[0];
	for (int i = 1; i < n; ++ i) {
		q[i-1] = c[i-1] / p[i-1];
		p[i] = a[i] - b[i-1] * q[i-1];
	}
	y[0] = y[0] / p[0];
	for (int i = 1; i < n; ++ i)
		y[i] = (y[i] - b[i-1] * y[i-1]) / p[i];
	for (int i = n-2; i >= 0 ; -- i)
		y[i] = y[i] - q[i] * y[i+1];
	return y;
}