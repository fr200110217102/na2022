#ifndef FUNCTION
#define FUNCTION
const double delta = 1e-8;
template <class type>
class Function{
public:
	virtual type operator ()(const type& x) const = 0;
	virtual type d(const type& x, const int& k = 1) const {
		if (k == 1) return ((*this)(x+delta) - (*this)(x-delta)) / (2*delta);
		if (k == 2) return ((*this)(x+delta) - 2*(*this)(x) + (*this)(x-delta)) / (delta*delta);
		throw 0;
	}
};
template <class type>
class Function_2D{
public:
	virtual type operator ()(const type& x, const type& y) const = 0;
	virtual type partial(const type& x, const type& y, const int& i, const int& j) const {
		if (i == 1 && j == 0) return ((*this)(x+delta, y) - (*this)(x-delta, y)) / (2*delta);
		if (i == 2 && j == 0) return ((*this)(x+delta, y) - 2*(*this)(x, y) + (*this)(x-delta, y))/ (delta*delta);
		if (i == 0 && j == 1) return ((*this)(x, y+delta) - (*this)(x, y-delta)) / (2*delta);
		if (i == 0 && j == 2) return ((*this)(x, y+delta) - 2*(*this)(x, y) + (*this)(x, y-delta))/ (delta*delta);
		throw 0;
	}
};
#endif