#include<bits/stdc++.h>
#include<Eigen/Core>
#include<Eigen/Dense>
#include<Eigen/LU>
#include<Eigen/QR>
#include<Eigen/Eigenvalues>
using namespace std;
using namespace Eigen;
int main(){
    vector<double> x = {0.0 ,0.5, 1.0, 1.5, 2.0, 2.5, 3.0,3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5,7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0};
    vector<double> y = {2.9 ,2.7 ,4.8 ,5.3 ,7.1 ,7.6, 7.7,7.6, 9.4, 9.0, 9.6, 10.0, 10.2, 9.7,8.3, 8.4, 9.0, 8.3, 6.6, 6.7, 4.1};
    MatrixXd a(21,3);
    VectorXd b(21);
    for(int i=0;i<21;++i){
        a(i,0)=1;
        a(i,1)=x[i];
        a(i,2)=x[i]*x[i];
        b(i)=y[i];
    }
    HouseholderQR<MatrixXd> qr(a);
    MatrixXd r=qr.matrixQR().triangularView<Upper>();
    MatrixXd q=qr.householderQ();
    cout<<"r = \n"<<r<<endl;
    MatrixXd r1=r.block<3,3>(0,0);
    VectorXd c=q.transpose()*b;
    VectorXd p=r1.lu().solve(c.head(3));
    cout<<"solution = \n"<<p<<endl;

	MatrixXd g(3,3);
	g<<21,105,717.5,105,717.5,5512.5,717.5,5512.5,45166.6;
	EigenSolver<MatrixXd> eig(g);
	cout<<"The Eigenvalues of G is:"<<endl<<eig.eigenvalues()<<endl;

	EigenSolver<MatrixXd> eig2(r1.transpose()*r1);
	cout<<"The Eigenvalues of R.T*R is:"<<endl<<eig2.eigenvalues()<<endl;
}