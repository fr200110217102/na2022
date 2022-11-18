#include<bits/stdc++.h>
#include"../HW1/function.h"
using namespace std;

class F : public Function<double> {
public:
	virtual double operator()(const double& x) const{
		return x*x*x*x*x*x*x*x-8*x*x*x*x*x*x*x+28*x*x*x*x*x*x-56*x*x*x*x*x+70*x*x*x*x-56*x*x*x+28*x*x-8*x+1;
	}
}f;

class G : public Function<double> {
public:
	virtual double operator()(const double& x) const{
		return (((((((x-8)*x+28)*x-56)*x+70)*x-56)*x+28)*x-8)*x+1;
	}
}g;

class H : public Function<double> {
public:
	virtual double operator()(const double& x) const{
		return (x-1)*(x-1)*(x-1)*(x-1)*(x-1)*(x-1)*(x-1)*(x-1);
	}
}h;

int main(){
	ofstream out("A.csv");
	out<<"x,f,g,h\n";
	for(int i=0;i<=100;++i){
		double x=0.99+i*0.0002;
		out<<x<<','<<f(x)<<','<<g(x)<<','<<h(x)<<'\n';
	}
}