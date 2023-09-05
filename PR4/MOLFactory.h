#include <bits/stdc++.h>
using namespace std;

#include <bits/stdc++.h>
#include "MOL.h"
using namespace std;
class MOLFactory {
public:
	using CreateCallback = unique_ptr<MOL>(*)();
private:
	using CallbackMap = map<string, CreateCallback>;
public:
	static MOLFactory& CreateFactory() {
		static MOLFactory object;
		return object;
	}
	bool Register(string ID, CreateCallback createFn) {
		return callbacks.insert({ID, createFn}).second;
	}
	bool UnRegister(string ID) {
		return callbacks.erase(ID) == 1;
	}
	unique_ptr<MOL> create(string ID) {
		auto it = callbacks.find(ID);
		if(it == callbacks.end()) {
			throw runtime_error("Unknown MOL ID. ");
		}
		return (it->second)();
	}
	private:
		MOLFactory() = default;
		MOLFactory(const MOLFactory &) = default;
		MOLFactory & operator = (const MOLFactory &) = default;
		~MOLFactory() = default;
	private:
		CallbackMap callbacks;
};