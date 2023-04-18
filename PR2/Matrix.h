#ifndef MATRIX
#define MATRIX
#include <bits/stdc++.h>
using namespace std;

const double eps = 1e-12;

template <class T>
class Matrix {
	protected:
		int n, m;
		T * a;
	public:
		Matrix() {
			a = NULL, n = m = 0;
		}
		void reserve(const int & n) {
			a = new T[n];
		}
		Matrix(int _n, int _m) : n(_n), m(_m) {
			reserve(n * m);
			memset(a, 0, sizeof(T) * n * m);
		}
		const T & at(const int & i) const {
			return a[i];
		}
		T & at(const int & i) {
			return a[i];
		}
		const T * operator [] (const int & i) const {
			return & a[i * m];
		}
		T * operator [] (const int & i) {
			return & a[i * m];
		}
		Matrix <T> (const Matrix <T> & p) {
			n = p.n, m = p.m;
			reserve(n * m);
			memcpy(a, p.a, sizeof(T) * n * m);
		}
		Matrix <T> (const initializer_list<initializer_list<T>>& args) {
			n = args.size(), m = args.begin()->size();
			reserve(n * m);
			for (int i = 0; i < n; ++ i) 
				for (int j = 0; j < m; ++ j)
					a[i*m+j] = *((*(args.begin() + i)).begin() + j);
		}
		Matrix <T> operator = (const Matrix <T> & p ) {
			if(a != NULL) delete a;
			n = p.n, m = p.m;
			reserve(n * m);
			memcpy(a, p[0], sizeof(T) * n * m);
			return *this;
		}
		~Matrix(){
			if(a != NULL) delete a;
		}
		inline int row() const {
			return n;
		}
		inline int col() const {
			return m;
		}
};

template <class T>
Matrix <T> operator + (const Matrix <T> & a, const Matrix <T> & b) {
	int n = a.row(), m = a.col();
	Matrix <T> c(n, m);
	for(int i = 0; i < n; ++ i)
		for(int j = 0; j < m; ++ j)
			c[i][j] = a[i][j] + b[i][j];
	return c;
}

template <class T>
Matrix <T> operator - (const Matrix <T> & a, const Matrix <T> & b) {
	int n = a.row(), m = a.col();
	Matrix <T> c(n, m);
	for(int i = 0; i < n; ++ i)
		for(int j = 0; j < m; ++ j)
			c[i][j] = a[i][j] - b[i][j];
	return c;
}

template <class T>
Matrix <T> operator * (const Matrix <T> & a, const T & t) {
	int n = a.row(), m = a.col();
	Matrix <T> c(n, m);
	for(int i = 0; i < n; ++ i)
		for(int j = 0; j < m; ++ j)
			c[i][j] = a[i][j] * t;
	return c;
}

template <class T>
Matrix <T> operator * (const T & t, const Matrix <T> & a) {
	return a * t;
}

template <class T>
Matrix <T> operator / (const Matrix <T> & a, const T & t) {
	return a * (1.0/t);
}

template <class T>
Matrix <T> operator - (const Matrix <T> & a) {
	return a * (-1.0);
}

template <class T>
Matrix <T> operator * (const Matrix <T> & a, const Matrix <T> & b) {
	int n = a.row(), m = a.col(), p = b.col();
	Matrix <T> c(n, p);
	for(int i = 0; i < n; ++ i)
		for(int j = 0; j < m; ++ j)
			for(int k = 0; k < p; ++ k)
				c[i][k] += a[i][j] * b[j][k];
	return c;
}

template <class T>
istream & operator >> (istream & in, Matrix <T> & a) {
	int n = a.row(), m = a.col();
	for(int i = 0; i < n; ++ i)
		for(int j = 0; j < m; ++ j)
			in >> a[i][j];
	return in;
}

template <class T>
ostream & operator << (ostream & out, const Matrix <T> & a) {
	int n = a.row(), m = a.col();
	if (n > 1 && m > 1) cout << "[";
	for(int i = 0; i < n; ++ i) {
		cout << "[";
		for(int j = 0; j < m-1; ++ j)
			out << a[i][j] << ", ";
		cout << a[i][m-1] << "]";
		if (i < n-1) cout << ",\n";
	}
	if (n > 1 && m > 1)cout << "]\n";
	return out;
}

template <class T>
Matrix <T> inverse(const Matrix <T> & a) {
	int n = a.row();
	Matrix<T> b(n, n * 2);
	for(int i = 0; i < n; ++ i)
		for(int j = 0; j < n; ++ j)
			b[i][j] = a[i][j];
	for(int i = 0; i < n; ++ i)
		b[i][i + n] = 1;
	for(int i = 0; i < n; ++ i) {
		int o = i;
		for(int j = i+1; j < n; ++ j)
			if(fabs(b[j][i]) > fabs(b[o][i])) o = j;
		if(o != i) {
			for(int j = 0; j < n*2; ++ j) swap(b[o][j], b[i][j]);
		}
		T tmp = b[i][i];
		for(int j = 0; j < n; ++ j) if(i ^ j) {
			T z = b[j][i] / b[i][i];
			for(int k = 0; k < n * 2; ++ k)
				b[j][k] -= b[i][k] * z;
		}
		for(int j = 0; j < n * 2;++ j)
			b[i][j] /= tmp;
	}
	Matrix<T> res(n, n);
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j)
				res[i][j] = b[i][j + n];
	return res;
}

template <class T>
Matrix <T> operator / (const Matrix <T> & a, const Matrix <T> & b) {
	return a * inverse(b);
}

template <class T>
Matrix <T> operator ~ (const Matrix <T> & a) {
	int n = a.row(), m = a.col();
	Matrix <T> b(m, n);
	for(int i = 0; i < n; ++ i)
		for(int j = 0; j < m; ++ j)
			b[j][i] = a[i][j];
	return b;
}

template <class T>
Matrix <T> Id (const int & n){
	Matrix <T> r(n, n);
	for(int i = 0; i < n; ++ i)
		r[i][i] = 1;
	return r;
}

template <class T> class Rowvec : public Matrix <T> {
	public:
		Rowvec() : Matrix <T> () {}
		Rowvec(int _n) : Matrix <T> (1, _n) {}
		Rowvec <T> (const Matrix <T> & p) : Matrix <T> (p) {}
		inline const T & operator [] (const int & i) const {
			return Matrix <T> :: at(i);
		}
		inline T & operator [] (const int & i) {
			return Matrix <T> :: at(i);
		}
};

template <class T> class Colvec : public Matrix <T> {
	public:
		Colvec() : Matrix <T> () {}
		Colvec(int _n) : Matrix <T> (_n, 1) {}
		Colvec <T> (const Matrix <T> &p) : Matrix <T> (p) {}
		Colvec <T> (const vector <T> &p) : Matrix <T> (p.size(), 1){
			int n = p.size();
			for (int i = 0; i < n; ++ i) Matrix <T> :: at(i) = p[i];
		}
		Colvec <T> (const initializer_list<T>& args) : Matrix <T>(args.size(), 1){
			int n = args.size();
			for (int i = 0; i < n; ++ i) Matrix <T> :: at(i) = *(args.begin() + i);
		}
		inline const T & operator [] (const int & i) const {
			return Matrix <T> :: at(i);
		}
		inline T & operator [] (const int & i) {
			return Matrix <T> :: at(i);
		}
};

template <class T>
Rowvec <T> operator + (const Rowvec <T> & a, const Rowvec <T> & b) {
	return (Matrix <T>) a + (Matrix <T> ) b;
}

template <class T>
Rowvec <T> operator - (const Rowvec <T> & a, const Rowvec <T> & b) {
	return (Matrix <T>) a - (Matrix <T> ) b;
}

template <class T>
Colvec <T> operator + (const Colvec <T> & a, const Colvec <T> & b) {
	return (Matrix <T>) a + (Matrix <T> ) b;
}

template <class T>
Colvec <T> operator - (const Colvec <T> & a, const Colvec <T> & b) {
	return (Matrix <T>) a - (Matrix <T> ) b;
}

template <class T>
Rowvec <T> operator * (const Rowvec <T> & a, const Matrix <T> & b) {
	return (Matrix <T>) a * b;
}

template <class T>
Colvec <T> operator * (const Matrix <T> & a, const Colvec <T> & b) {
	return a * (Matrix <T>) b;
}

template <class T>
Colvec <T> operator * (const T & c, const Colvec <T> & b) {
	return c * (Matrix <T>) b;
}

template <class T>
Colvec <T> operator * (const Colvec <T> & b, const T & c) {
	return c * b;
}

template <class T>
inline T operator * (const Rowvec <T> & a, const Colvec <T> & b) {
	return (Matrix<T>(a) * Matrix<T>(b)).at(0);
}

template <class T>
inline Rowvec <T> operator ~ (const Colvec <T> & a) {
	return (Rowvec<T>)(~(Matrix <T>(a)));
}

template <class T>
inline Colvec <T> operator ~ (const Rowvec <T> & a) {
	return (Colvec<T>)(~(Matrix <T>(a)));
}

template <class T>
inline Colvec <T> split(const Colvec <T> & a, int l, int r) {
	Colvec <T> b(r - l);
	for(int i = l; i < r; ++ i)
		b[i-l] = a[i];
	return b;
}

template <class T>
inline Colvec <T> split(const Matrix <T> & a, int p, int l, int r) {
	Colvec <T> b(r - l);
	for(int i = l; i < r; ++ i)
		b[i-l] = a[i][p];
	return b;
}

template <class T>
inline Matrix <T> split(const Matrix <T> & a, int u, int d, int l, int r) {
	Matrix <T> b(d - u, r - l);
	for(int i = u; i < d; ++ i)
		for(int j = l; j < r; ++ j)
			b[i-u][j-l] = a[i][j];
	return b;
}

template <class T>
void update(Colvec <T> & a, const Colvec <T> & b, int l, int r) {
	for(int i = l; i < r; ++ i)
		a[i] = b[i-l];
}

template <class T>
void update(Matrix <T> & a, const Matrix <T> & b, int u, int d, int l, int r) {
	for(int i = u; i < d; ++ i)
		for(int j = l; j < r; ++ j)
			a[i][j] = b[i-u][j-l];
}

template <class T>
Colvec <T> e(const int & n, const int & i) {
	Colvec <T> x(n);
	x[i] = 1;
	return x;
}

template <class T>
Colvec <T> FSP(const Matrix <T> & a, const Colvec <T> & b, const int & op = 1) {
	int n = a.row();
	Colvec <T> x = b;
	for(int i = 0; i < n; ++ i){
		for(int j = 0; j < i; ++ j) x[i] -= a[i][j] * x[j];
		if(op) x[i] /= a[i][i];
	}
	return x;
}

template <class T>
Colvec <T> BSP(const Matrix <T> & a, const Colvec <T> & b, const int & op = 1) {
	int n = a.row();
	Colvec <T> x = b;
	for(int i = n-1; ~i; -- i){
		for(int j = n-1; j > i; -- j) x[i] -= a[i][j] * x[j];
		if(op) x[i] /= a[i][i];
	}
	return x;
}

template <class T>
Matrix <T> Gauss_LU(const Matrix <T> & _a) {
	Matrix <T> a = _a;
	int n = _a.row();
	for(int k = 0; k < n; ++ k)
		for(int i = k+1; i < n; ++ i) {
			a[i][k] /= a[k][k];
			for(int j = k+1; j < n; ++ j)
				a[i][j] -= a[i][k] * a[k][j];
		}
	return a;
}

template <class T>
Colvec <T> Gauss_Solve(const Matrix <T> & a, const Colvec <T> & b) {
	Matrix <T> u = Gauss_LU(a);
	return BSP(u, FSP(u, b, 0));
}

template <class T>
pair <Matrix <T>, vector<int> > Gauss_Improved(const Matrix <T> & _a) {
	Matrix <T> a = _a;
	int n = _a.row();
	vector <int> p(n);
	for(int i = 0; i < n; ++ i) p[i] = i;
	for(int k = 0; k < n; ++ k){
		int o = k;
		for(int i = k+1; i < n; ++ i)
			if(fabs(a[i][k]) > fabs(a[o][k])) o = i;
		if(o != k) {
			for(int j = 0; j < n; ++ j) swap(a[o][j], a[k][j]);
			swap(p[k], p[o]);
		}
		vector<int> cols;
		for(int j = k; j < n; ++ j) if (a[k][j]) cols.push_back(j);
		for(int i = k+1; i < n; ++ i) {
			if (!a[i][k]) continue;
			a[i][k] /= a[k][k];
			for(int j = 1; j < cols.size(); ++ j)
				a[i][cols[j]] -= a[i][k] * a[k][cols[j]];
		}
	}
	return make_pair(a, p);
}

template <class T>
Colvec <T> Gauss_Improved_Solve(const Matrix <T> & _a, const Colvec <T> & _b){
	pair <Matrix <T>, vector <int> > gs = Gauss_Improved(_a);
	Colvec <T> b = _b;
	int n = _a.row();
	for(int i = 0; i < n; ++ i) b[i] = _b[gs.second[i]];
	return BSP(gs.first, FSP(gs.first, b, 0));
}

template <class T>
Matrix <T> Cholesky_Improved(const Matrix <T> & a) {
	int n = a.row();
	Matrix <T> b = a;
	for(int k = 0; k < n; ++ k) {
		for(int j = 0; j < k; ++ j) b[k][k] -= b[k][j] * b[k][j] * b[j][j];
		for(int i = k+1; i < n; ++ i){
			for(int j = 0; j < k; ++ j)
				b[i][k] -= b[i][j] * b[j][j] * b[k][j];
			b[i][k] /= b[k][k];
		}
	}
	return b;
}

template <class T>
Colvec <T> Cholesky_Improved_Solve(const Matrix <T> & a, const Colvec <T> & b) {
	int n = a.row();
	Matrix <T> l = Cholesky_Improved(a);
	Colvec <T> x = FSP(l, b, 0);
	for(int i = 0; i < n; ++ i) x[i] /= l[i][i];
	return BSP(~l, x, 0);
}

template <class T>
double vert_1(const Colvec <T> & a) {
	int n = a.row();
	double res = 0;
	for(int i = 0; i < n; ++ i) res += abs(a[i]);
	return res;
}

template <class T>
double vert_2_2(const Colvec <T> & a) {
	int n = a.row();
	double res = 0;
	for(int i = 0; i < n; ++ i) res += abs(a[i] * a[i]);
	return res;
}

template <class T>
double vert_2(const Colvec <T> & a) {
	return sqrt(vert_2_2(a));
}

template <class T>
double vert_inf(const Colvec <T> & a) {
	int n = a.row();
	double res = 0;
	for(int i = 0; i < n; ++ i) res = max(res, abs(a[i]));
	return res;
}

template <class T>
double Vert_1(const Matrix <T> & a) {
	int n = a.row(), m = a.col();
	T res = 0;
	for(int j = 0; j < m; ++ j) {
		T s = 0;
		for(int i = 0; i < n; ++ i)
			s += abs(a[i][j]);
		res = max(res, s);
	}
	return res;
}

template <class T>
double Vert_inf(const Matrix <T> & a) {
	return Vert_1(~a);
}

template <class T>
inline int sgn(const T & x) {
	return x >= 0 ? 1 : -1;
}

template <class T>
Colvec <T> sgn(const Colvec <T> & a) {
	int n = a.row();
	Colvec <T> s(n);
	for(int i = 0; i < n; ++ i) s[i] = sgn(a[i]);
	return s;
}

template <class T>
double Kappa(const Matrix <T> & a) {
	int n = a.row();
	Matrix <T> u = Gauss_LU(a);
	Colvec <T> x = e<T>(n, 0), w(n), v(n), z(n);
	while(1) {
		w = BSP(~u, FSP(~u, x), 0);
		v = sgn(w);
		z = BSP(u, FSP(u, v, 0));
		T d = vert_inf(z);
		if(d <= (~z)*x) return Vert_inf(a) * vert_1(w);
		int p = 0;
		for(int i = 1; i < n; ++ i) if(fabs(z[i]) > fabs(z[p])) p = i;
		x = e<T>(n, p);
	}
}

template <class T>
pair <Colvec <T>, T> LS(const Matrix <T> & a, const Colvec <T> & b) {
	Colvec <T> x = Cholesky_Improved_Solve((~a)*a, (~a)*b);
	T r = vert_2(a * x - b);
	return make_pair(x, r);
}

template <class T>
pair <Colvec <T>, T> house(const Colvec <T> & x) {
	int n = x.row();
	if (fabs(vert_inf(x)) < eps) return make_pair(e<T>(n, 0), 0);
	T S = 0;
	Colvec <T> v = x / vert_inf(x);
	for(int i = 1; i < n; ++ i) S += v[i] * v[i];
	if (S < eps) return make_pair(v, 0);
	T A = sqrt(S + v[0] * v[0]);
	if (v[0] <= 0) v[0] -= A;
	else v[0] = -S/(v[0] + A);
	T B = 2 * v[0] * v[0] / (S + v[0] * v[0]);
	return make_pair(v / v[0], B);
}

template <class T>
Matrix <T> house_trans(const Matrix <T> & a, const pair <Colvec <T>, T> & vb, bool type = 0) {
	if (!type) {
		Colvec <T> v = vb.first;
		v[0] = 1;
		return a - v * ~(~a * v * vb.second);
	} else {
		Colvec <T> v = vb.first;
		v[0] = 1;
		return a - a * (vb.second * v) * ~v;
	}
}

template <class T>
pair <Matrix <T>, Colvec<T> > QR(const Matrix <T> & a) {
	Matrix <T> b = a;
	int m = a.row(), n = a.col();
	Colvec <T> d(n);
	pair <Colvec <T>, T> h;
	for(int k = 0; k < n; ++ k){
		h = house(split(b, k, k, m));
		update(b, house_trans(split(b, k, m, k, n), h), k, m, k, n);
		d[k] = h.second;
		for(int i = k+1; i < m; ++ i) b[i][k] = h.first[i-k];
	}
	return make_pair(b, d);
}

template <class T>
pair <Colvec <T>, T> QR_LS(const Matrix <T> & a, const Colvec <T> & b) {
	int m = a.row(), n = a.col();
	pair <Matrix <T>, Colvec <T> > qr = QR(a);
	Colvec <T> c = b;
	for(int k = 0; k < n; ++ k)
		update(c, (Colvec <T>)house_trans(split(c, k, m), make_pair(split(qr.first, k, k, m), qr.second[k])), k, m);
	Colvec <T> x = BSP(split(qr.first, 0, n, 0, n), split(c, 0, n));
	T res = vert_2(split(c, n, m));
	return make_pair(x, res);
}

template <class T>
Colvec <T> Jacobi(const Matrix <T> & a, const Colvec <T> & b, const T & eps = 1e-12) {
	int n = a.row(), times = 0;
	Colvec <T> now = b, last(n);
	do{
		++ times;
		last = now;
		for(int i = 0; i < n; ++ i) {
			now[i] = 0;
			for(int j = 0; j < n; ++ j) if(j != i)
				now[i] += a[i][j] * last[j];
			now[i] -= b[i];
			now[i] /= -a[i][i];
		}
	}while(vert_2(last - now) >= eps);
	// cout << "times = " << times << endl;
	return now;
}

template <class T>
Colvec <T> Gauss_Seidel(const Matrix <T> & a, const Colvec <T> & b, const T & eps = 1e-12) {
	int n = a.row(), times = 0;
	Colvec <T> now = b, last(n);
	do{
		++ times;
		last = now;
		for(int i = 0; i < n; ++ i) {
			now[i] = 0;
			for(int j = 0; j < i; ++ j)
				now[i] += a[i][j] * now[j];
			for(int j = i+1; j < n; ++ j)
				now[i] += a[i][j] * last[j];
			now[i] -= b[i];
			now[i] /= -a[i][i];
		}
	}while(vert_2(last - now) >= eps);
	cout << "times = " << times << endl;
	return now;
}

template <class T>
Colvec <T> SOR(const Matrix <T> & a, const Colvec <T> & b, const T & w, const T & eps = 1e-12) {
	int n = a.row(), times = 0;
	Colvec <T> now = b, last(n);
	do{
		++ times;
		last = now;
		for(int i = 0; i < n; ++ i) {
			now[i] = 0;
			for(int j = 0; j < i; ++ j)
				now[i] += a[i][j] * now[j];
			for(int j = i+1; j < n; ++ j)
				now[i] += a[i][j] * last[j];
			now[i] -= b[i];
			now[i] *= -(w/a[i][i]);
			now[i] += (1-w) * last[i];
		}
	}while(vert_2(last - now) >= eps);
	// cout << "times = " << times << endl;
	return now;
}

template <class T>
Colvec <T> CG(const Matrix <T> & a, const Colvec <T> & b, const T & eps = 1e-12, const int & maxtime = 20000) {
	int n = a.row(), times = 0;
	Colvec <T> x(n), r(n), p(n), w(n);
	T alpha, beta, rho, lastrho;
	r = b - a * x;
	rho = (~r) * r;
	while(vert_2(r) >= eps && times < maxtime) {
		++ times;
		if(times == 1) p = r;
		else{
			beta = rho / lastrho;
			p = r + beta * p;
		}
		w = a * p;
		alpha = rho / ((~p) * w);
		x = x + alpha * p;
		r = r - alpha * w;
		lastrho = rho;
		rho = (~r) * r;
	}
	// cout << "times = " << times << endl;
	return x;
}

template <class T>
Matrix <T> Hessenberg(Matrix <T> & a) {
	int n = a.row();
	Matrix <T> Q = Id <T>(n);
	pair <Colvec <T>, T> vb;
	for (int k = 0; k < n-2; ++ k) {
		vb = house(split(a, k, k+1, n));
		update(a, house_trans(split(a, k+1, n, k, n), vb, 0), k+1, n, k, n);
		update(a, house_trans(split(a, 0, n, k+1, n), vb, 1), 0, n, k+1, n);
		update(Q, house_trans(split(Q, 0, n, k+1, n), vb, 1), 0, n, k+1, n);
	}
	return Q;
}

template <class T>
vector < pair <Colvec <T>, T> > QR_iteration(Matrix <T> & h) {
	vector <pair <Colvec <T>, T> > res(0);
	pair <Colvec <T>, T> vb;
	Colvec <T> tmp(3);
	int n = h.row();
	if (n == 2) {
		T mu = h[1][1];
		h[0][0] -= mu, h[1][1] = 0;
		T x = mu;
		T y = h[1][0];
		tmp = Colvec <T> (2);
		tmp[0] = x, tmp[1] = y;
		res.push_back(vb = house(tmp));
		h = house_trans(h, vb, 0);
		h = house_trans(h, vb, 1);
		h[0][0] += mu, h[1][1] += mu;
		return res;
	}
	T s = h[n-2][n-2] + h[n-1][n-1];
	T t = h[n-2][n-2] * h[n-1][n-1] - h[n-2][n-1] * h[n-1][n-2];
	T x = h[0][0] * h[0][0] + h[0][1] * h[1][0] - s * h[0][0] + t;
	T y = h[1][0] * (h[0][0] + h[1][1] - s);
	T z = h[1][0] * h[2][1];
	for (int k = 0; k < n-2; ++ k) {
		tmp[0] = x, tmp[1] = y, tmp[2] = z;
		res.push_back(vb = house(tmp));
		int q = max(0, k-1);
		update(h, house_trans(split(h, k, k+3, q, n), vb, 0), k, k+3, q, n);
		int r = min(k+4, n);
		update(h, house_trans(split(h, 0, r, k, k+3), vb, 1), 0, r, k, k+3);
		x = h[k+1][k];
		y = h[k+2][k];
		if (k < n-3) z = h[k+3][k];
	}
	tmp = Colvec <T> (2);
	tmp[0] = x, tmp[1] = y;
	res.push_back(vb = house(tmp));
	update(h, house_trans(split(h, n-2, n, n-3, n), vb, 0), n-2, n, n-3, n);
	update(h, house_trans(split(h, 0, n, n-2, n), vb, 1), 0, n, n-2, n);
	return res;
}

template <class T>
inline complex <T> eigen_2(const T& a, const T& b, const T& c, const T& d) {
	T p = a + d, q = a * d - b * c, delta = p * p - 4 * q;
	if (delta < 0) return complex <T> (p/2, sqrt(-delta)/2);
	else return complex <T> ((p + sqrt(delta)) / 2, 0);
}

template <class T>
Matrix <T> Schur(Matrix <T> & a, const T & eps = 1e-12) {
	Matrix <T> h22;
	int n = a.row(), l, r, m, n0;
	vector < pair < Colvec <T>, T> > trs;
	Matrix <T> Q = Hessenberg(a);
	int times = 0;
	while (1) {
		for (int i = 1; i < n; ++ i)
			if (fabs(a[i][i-1]) <= (fabs(a[i][i]) + fabs(a[i-1][i-1])) * eps)
				a[i][i-1] = 0;
		r = n-1;
		while (r >= 1) {
			if (fabs(a[r][r-1]) < eps) -- r;
			else if ((r == 1 || (r > 1 && fabs(a[r-1][r-2]) < eps)) && fabs(eigen_2(a[r-1][r-1], a[r-1][r], a[r][r-1], a[r][r]).imag()) > eps) r -= 2;
			else break;
		}
		if (r <= 0) break;
		l = r, ++ r, m = n-r;
		while (l >= 1) {
			if (fabs(a[l][l-1]) > eps) -- l;
			else break;
		}
		n0 = r-l;
		h22 = split(a, l, r, l, r);
		trs = QR_iteration(h22);
		
		for (int k = 0; k < n0-2; ++ k) {
			if (l) update(a, house_trans(split(a, 0, l, l+k, l+k+3), trs[k], 1), 0, l, l+k, l+k+3);
			if (m) update(a, house_trans(split(a, l+k, l+k+3, r, n), trs[k], 0), l+k, l+k+3, r, n);
			update(Q, house_trans(split(Q, 0, n, l+k, l+k+3), trs[k], 1), 0, n, l+k, l+k+3);
		}
		if (l) update(a, house_trans(split(a, 0, l, r-2, r), trs[n0-2], 1), 0, l, r-2, r);
		if (m) update(a, house_trans(split(a, r-2, r, r, n), trs[n0-2], 0), r-2, r, r, n);
		update(Q, house_trans(split(Q, 0, n, r-2, r), trs[n0-2], 1), 0, n, r-2, r);
		
		update(a, h22, l, r, l, r);
		++ times;
	}
	// cout << "times = " << times << endl;
	return Q;
}

template <class T>
vector < complex <T> > QR_eigens(Matrix <T> & a) {
	Schur(a);
	int n = a.row(), r = 0;
	vector <complex <T> > eigens;
	while (r < n-1) {
		if (fabs(a[r+1][r]) < 1e-12) eigens.push_back(a[r][r]), ++r;
		else {
			complex <double> lambda = eigen_2(a[r][r], a[r][r+1], a[r+1][r], a[r+1][r+1]);
			eigens.push_back(lambda), eigens.push_back(conj(lambda));
			r += 2;
		}
	}
	if (r == n-1) eigens.push_back(a[r][r]);
	return eigens;
}

template <class T>
Matrix <T> Symmetric_Hessenberg(Matrix <T> & a) {
	int n = a.row();
	Matrix <T> Q = Id <T> (n);
	pair <Colvec <T>, T> vb;
	Colvec <T> u, v, w;
	T bt;
	for (int k = 0; k < n-2; ++ k) {
		vb = house(split(a, k, k+1, n));
		v = vb.first, bt = vb.second;
		u = bt * (split(a, k+1, n, k+1, n) * v);
		w = u - (bt * (~u * v) / 2.0) * v;
		a[k][k+1] = a[k+1][k] = vert_2(split(a, k, k+1, n));
		Matrix <T> tmp = split(a, k+1, n, k+1, n) - v * ~w - w * ~v;
		update(a, tmp, k+1, n, k+1, n);
		update(Q, house_trans(split(Q, 0, n, k+1, n), vb, 1), 0, n, k+1, n);
	}
	return Q;
}

template <class T>
vector <pair<T, T> > Symmetric_QR_iteration(Matrix <T> & t) {
	int n = t.row();
	T d = (t[n-2][n-2] - t[n-1][n-1]) / 2;
	T mu = t[n-1][n-1] - t[n-1][n-2] * t[n-1][n-2] / (d + sgn(d) * sqrt(d * d + t[n-1][n-2] * t[n-1][n-2]));
	T x = t[0][0] - mu, z = t[1][0];
	Matrix <T> t1, t2;
	vector <pair<T, T> > res(0);
	for (int k = 0; k < n-1; ++ k) {
		T c = x / sqrt(x * x + z * z), s = z / sqrt(x * x + z * z);
		int l = max(0, k-1), r = min(k+4, n);
		t1 = split(t, k,   k+1, l, r);
		t2 = split(t, k+1, k+2, l, r);
		update(t, t1 * c  + t2 * s, k,   k+1, l, r);
		update(t, t1 * -s + t2 * c, k+1, k+2, l, r); 
		t1 = split(t, l, r, k,   k+1);
		t2 = split(t, l, r, k+1, k+2);
		update(t, t1 * c  + t2 * s, l, r, k,   k+1);
		update(t, t1 * -s + t2 * c, l, r, k+1, k+2); 
		res.push_back(make_pair(c, s));
		if (k < n-2) x = t[k+1][k], z = t[k+2][k];
	}
	return res;
}

template <class T>
Matrix <T> Symmetric_Diagnize(Matrix <T> & a, const T eps = 1e-12) {
	Matrix <T> t22;
	int n = a.row(), l, r, m, n0;
	vector <pair<T, T> > trs;
	Matrix <T> Q = Symmetric_Hessenberg(a);
	Matrix <T> t1, t2;
	int times = 0;
	while (1) {
		for (int i = 1; i < n; ++ i)
			if (fabs(a[i][i-1]) <= (fabs(a[i][i]) + fabs(a[i-1][i-1])) * eps)
				a[i][i-1] = 0;
		r = n-1;
		while (r >= 1) {
			if (fabs(a[r][r-1]) < eps) -- r;
			else break;
		}
		if (r <= 0) break;
		l = r, ++ r, m = n-r;
		while (l >= 1) {
			if (fabs(a[l][l-1]) > eps) -- l;
			else break;
		}
		n0 = r-l;
		t22 = split(a, l, r, l, r);
		trs = Symmetric_QR_iteration(t22);
		
		for (int k = 0; k < n0-1; ++ k) {
			t1 = split(Q, 0, n, l+k  , l+k+1);
			t2 = split(Q, 0, n, l+k+1, l+k+2);
			T c = trs[k].first, s = trs[k].second;
			update(Q, t1 *  c + t2 * s, 0, n, l+k,   l+k+1);
			update(Q, t1 * -s + t2 * c, 0, n, l+k+1, l+k+2);
		}
		
		update(a, t22, l, r, l, r);
		++ times;
	}
	// cout << "times = " << times << endl;
	return Q;
}

template <class T> 
Matrix <T> Jacobi_Diagnize(Matrix <T> & a, const T eps1 = 1e-10, const T eps2 = 1e-12) {
	int n = a.row();
	Matrix <T> Q = Id <T>(n), tmp(n, n);
	int times = 0;
	while (1) {
		++ times;
		for (int p = 0; p < n-1; ++ p)
			for (int q = p+1; q < n; ++ q) if (fabs(a[p][q]) > eps2) {
				T r = (a[q][q] - a[p][p]) / (2 * a[p][q]);
				T t = 1.0 / (r >= 0 ? r + sqrt(1 + r * r) : r - sqrt(1 + r * r));
				T c = 1.0 / sqrt(1 + t * t), s = t / sqrt(1 + t * t);
				for (int i = 0; i < n; ++ i) if (i != p && i != q) {
					tmp[i][p] = c * a[i][p] - s * a[i][q];
					tmp[i][q] = c * a[i][q] + s * a[i][p];
				}
				tmp[p][p] = a[p][p] - t * a[p][q];
				tmp[q][q] = a[q][q] + t * a[p][q];
				tmp[p][q] = tmp[q][p] = 0;
				for (int i = 0; i < n; ++ i)
					a[i][p] = a[p][i] = tmp[i][p], a[i][q] = a[q][i] = tmp[i][q];
				for (int i = 0; i < n; ++ i) {
					tmp[i][p] =  c * Q[i][p] - s * Q[i][q];
					tmp[i][q] =  c * Q[i][q] + s * Q[i][p];
				}
				for (int i = 0; i < n; ++ i)
					Q[i][p] = tmp[i][p], Q[i][q] = tmp[i][q];
			}
		T sum = 0;
		for (int i = 0; i < n; ++ i)
			for (int j = 0; j < n; ++ j) if(i != j)
				sum += a[i][j] * a[i][j];
		if (sum < eps1 * n) break;
	}
	// cout << "times = " << times << endl;
	return Q;
}

template <class T>
T V(Colvec <T> & a, Colvec <T> & b, const T x) {
	int n = a.row();
	T now = a[0] - x;
	int res = now < 0;
	for (int i = 1; i < n; ++ i) {
		now = a[i] - x - b[i-1] * b[i-1] / now;
		res += now < 0;
	}
	return res;
}

template <class T>
T mth_eigen(Colvec <T> & a, Colvec <T> & b, int m, const T eps = 1e-12) {
	int n = a.row();
	T l, r = 0, mid;
	r = max(fabs(a[0]) + fabs(b[0]), fabs(a[n-1]) + fabs(b[n-2]));
	for (int i = 1; i < n-1; ++ i) r = max(r, fabs(a[i]) + fabs(b[i-1]) + fabs(b[i]));
	l = -r;
	while (r - l > eps) {
		mid = (l + r) / 2;
		if (V(a, b, mid) >= m) r = mid;
		else l = mid;
	}
	return mid;
}

template <class T>
T calc_prod(Colvec <T> & a, Colvec <T> & b, const Colvec <T> & x) {
	int n = a.row();
	T res = 0, w;
	w = a[0] * x[0] + b[0] * x[1];
	res += w * w;
	for (int i = 1; i < n-1; ++ i) {
		w = a[i] * x[i] + b[i-1] * x[i-1] + b[i] * x[i+1];
		res += w * w;
	}
	w = a[n-1] * x[n-1] + b[n-2] * x[n-2];
	res += w * w;
	return sqrt(res);
}

template <class T>
Colvec <T> eigen_vec(Colvec <T> & a, Colvec <T> & b, const T lmd, const T eps = 1e-8) {
	int n = a.row();
	Colvec <T> r(n), c(n), d(n-1), x(n), y(n);
	for (int i = 0; i < n; ++ i) r[i] = a[i] - lmd;
	c[0] = r[0], d[0] = b[0] / r[0];
	for (int k = 1; k < n-1; ++ k) {
		c[k] = r[k] - b[k-1] * d[k-1];
		d[k] = b[k] / c[k];
	}
	c[n-1] = r[n-1] - b[n-2] * d[n-2];
	x = e <T> (n, 0), y = x;
	while (calc_prod(r, b, y) > eps) {
		x = y;
		for (int k = 0; k < n; ++ k) {
			if (k > 0) x[k] -= x[k-1] * b[k-1];
			x[k] /= c[k];
		}
		for (int k = n-2; k >= 0; -- k)
			x[k] -= x[k+1] * d[k];
		y = x / vert_inf(x);
	}
	return y;
}
#endif