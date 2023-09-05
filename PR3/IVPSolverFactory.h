#include <bits/stdc++.h>
#include "IVPSolver.h"
using namespace std;
class IVPSolverFactory {
public:
	using CreateCallback = unique_ptr<IVPSolver>(*)();
private:
	using CallbackMap = map<string, CreateCallback>;
public:
	static IVPSolverFactory& CreateFactory() {
		static IVPSolverFactory object;
		return object;
	}
	bool Register(string ID, CreateCallback createFn) {
		return callbacks.insert({ID, createFn}).second;
	}
	bool UnRegister(string ID) {
		return callbacks.erase(ID) == 1;
	}
	unique_ptr<IVPSolver> create(string ID) {
		auto it = callbacks.find(ID);
		if(it == callbacks.end()) {
			throw runtime_error("Unknown IVPSolver ID. ");
		}
		return (it->second)();
	}
	private:
		IVPSolverFactory() = default;
		IVPSolverFactory(const IVPSolverFactory &) = default;
		IVPSolverFactory & operator = (const IVPSolverFactory &) = default;
		~IVPSolverFactory() = default;
	private:
		CallbackMap callbacks;
};