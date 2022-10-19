#ifndef FUNCTION
#define FUNCTION
const double delta = 1e-12;
template <class type>
class Function{
public:
	virtual type operator ()(const type& x) const = 0;
	virtual type d(const type& x) const {
		return ((*this)(x+delta) - (*this)(x-delta)) / (2*delta);
	}
};
#endif