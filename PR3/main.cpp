#include <bits/stdc++.h>
#include "Matrix.h"
#include "Function.h"
#include "LMM.h"
#include "RK.h"
#include "Emb_RK.h"
#include "IVPSolverFactory.h"
using namespace std;

class F : public Function_nd<double>{
private:
	double mu;
public:
	F(double mu):mu(mu){n=6;}
	Colvec<double> operator()(const Colvec<double>& u) const {
		if (u.row() != 6) throw "Error: Wrong Dimension of Vector!";
		return {
			u[3],
			u[4],
			u[5],
			2*u[4]+u[0]
			-mu*(u[0]+mu-1)		/pow(u[1]*u[1]+u[2]*u[2]+(u[0]+mu-1)*(u[0]+mu-1),1.5)
			-(1-mu)*(u[0]+mu)	/pow(u[1]*u[1]+u[2]*u[2]+(u[0]+mu)*(u[0]+mu)	,1.5),
			-2*u[3]+u[1]
			-mu*u[1]			/pow(u[1]*u[1]+u[2]*u[2]+(u[0]+mu-1)*(u[0]+mu-1),1.5)
			-(1-mu)*u[1]		/pow(u[1]*u[1]+u[2]*u[2]+(u[0]+mu)*(u[0]+mu)	,1.5),
			-mu*u[2]			/pow(u[1]*u[1]+u[2]*u[2]+(u[0]+mu-1)*(u[0]+mu-1),1.5)
			-(1-mu)*u[2]		/pow(u[1]*u[1]+u[2]*u[2]+(u[0]+mu)*(u[0]+mu)	,1.5)
		};
	}
};

double T[2] = {17.06521656015796, 19.140540691377};
Colvec<double> u0[2] = {
	{0.994, 0, 0, 0, -2.0015851063790825224, 0},
	{0.879779227778, 0, 0, 0, -0.379677780949, 0}
};

int main(int argc, char** argv) {
	IVPSolverFactory& Factory = IVPSolverFactory::CreateFactory();
	Factory.Register("Adam-Bashforth", []() -> unique_ptr<IVPSolver>{return make_unique<Adam_Bashforth>();});
	Factory.Register("Adam-Moulton", []() -> unique_ptr<IVPSolver>{return make_unique<Adam_Moulton>();});
	Factory.Register("BDF", []() -> unique_ptr<IVPSolver>{return make_unique<BDF>();});
	Factory.Register("Classical-RK", []() -> unique_ptr<IVPSolver>{return make_unique<Classical_4th_RK>();});
	Factory.Register("ESDIRK", []() -> unique_ptr<IVPSolver>{return make_unique<ESDIRK>();});
	Factory.Register("Gauss-Legendre-1", []() -> unique_ptr<IVPSolver>{return make_unique<Gauss>(1);});
	Factory.Register("Gauss-Legendre-2", []() -> unique_ptr<IVPSolver>{return make_unique<Gauss>(2);});
	Factory.Register("Gauss-Legendre-3", []() -> unique_ptr<IVPSolver>{return make_unique<Gauss>(3);});
	Factory.Register("Fehlberg", []() -> unique_ptr<IVPSolver>{return make_unique<Fehlberg_Emb_RK>();});
	Factory.Register("Dormand-Prince", []() -> unique_ptr<IVPSolver>{return make_unique<Dormand_Prince_Emb_RK>();});

	string name = argv[1];
	if (name == "Gauss-Legendre") name = name + '-' + argv[2];
	int s = atoi(argv[2]), n = atoi(argv[3]), type = atoi(argv[4]);

	F f(0.012277471);
	unique_ptr<IVPSolver> Solver = Factory.create(name);
	vector<Colvec<double>> u = Solver->Solve(f, u0[type], T[type], n, s);

	if (argc > 4) {
		ofstream out(argv[5], ios::out);
		out << "x,y\n";
		for (int i = 0; i < u.size(); ++ i) out << u[i][0] << ',' << u[i][1] << '\n';
		out.flush();
		out.close();
	}
}