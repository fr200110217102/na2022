#include<bits/stdc++.h>
#include "DLS.h"
using namespace std;
int main(){
	vector<double> x({0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10});
	vector<double> y({2.9,2.7,4.8,5.3,7.1,7.6,7.7,7.6,9.4,9.0,9.6,10.0,10.2,9.7,8.3,8.4,9.0,8.3,6.6,6.7,4.1});
	DLS<double> dlsn(x,y,2);
	DLS<double> dlsq(x,y,2,"QR");
	cout<<dlsn.getpoly()<<endl;
	cout<<dlsq.getpoly()<<endl;
}