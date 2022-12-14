#include <bits/stdc++.h>
#include "../HW2/Polynomial.h"
#include "../Matrix.h"
using namespace std;

template<class type>
class DLS{
private:
	Polynomial <type> p;
public:
	DLS(vector<type>& x,vector<type>& y,int m,string method="Normal"){
		if(x.size()!=y.size())throw "Invalid Initial!";
		if(method!="Normal"&&method!="QR")throw "Invalid Method!";
		int n=x.size();
		p.resize(m+1);
		if(method=="Normal"){
			vector<vector<type>> px(2*m+1),py(m+1);
			vector<type> pxs(2*m+1),pys(m+1);
			Matrix<type> G(m+1,m+1);
			Colvec<type> c(m+1),sol;
			px[0].resize(n);
			for(int i=0;i<n;++i)px[0][i]=1;
			pxs[0]=n;
			for(int i=1;i<=2*m;++i){
				px[i].resize(n);
				for(int j=0;j<n;++j){
					px[i][j]=px[i-1][j]*x[j];
					pxs[i]+=px[i][j];
				}
			}
			py[0].resize(n);
			for(int i=0;i<n;++i)py[0][i]=y[i],pys[0]+=y[i];
			for(int i=1;i<=m;++i){
				py[i].resize(n);
				for(int j=0;j<n;++j){
					py[i][j]=py[i-1][j]*x[j];
					pys[i]+=py[i][j];
				}
			}
			for(int i=0;i<=m;++i)
				for(int j=0;j<=m;++j)
					G[i][j]=pxs[i+j];
			for(int i=0;i<=m;++i)c[i]=pys[i];
			sol=Gauss_Improved_Solve(G,c);
			for(int i=0;i<=m;++i)p[i]=sol[i];
			/*答题用，使用时注释掉*/
			cout<<G<<endl;
			Jacobi_Diagnize(G);
			vector<type> eig(m+1);
			for(int i=0;i<=m;++i)eig[i]=fabs(G[i][i]),cout<<eig[i]<<' ';
			cout<<*max_element(eig.begin(),eig.end())<<endl<<endl;
			/**/
		}
		else if(method=="QR"){
			p.resize(m+1);
			Matrix<type> A(n,m+1);
			Colvec<type> b(n);
			for(int i=0;i<n;++i){
				A[i][0]=1;
				for(int j=1;j<=m;++j)A[i][j]=A[i][j-1]*x[i];
			}
			for(int i=0;i<n;++i)b[i]=y[i];
			pair<Matrix<type>,Colvec<type>> qr=QR(A);
			/*答题用，使用时注释掉*/
			Matrix<type> R(m+1,m+1);
			for(int i=0;i<=m;++i){
				for(int j=0;j<i;++j)R[i][j]=0;
				for(int j=i;j<=m;++j)R[i][j]=qr.first[i][j];
			}
			cout<<R<<endl;
			R=~R*R;
			Jacobi_Diagnize(R);
			vector<type> eig(m+1);
			for(int i=0;i<=m;++i)eig[i]=fabs(R[i][i]),cout<<eig[i]<<' ';
			cout<<sqrt(*max_element(eig.begin(),eig.end()))<<endl;
			/**/
			pair<Colvec<type>,type> sol=QR_LS(A,b);
			for(int i=0;i<=m;++i)p[i]=sol.first[i];
		}
	}
	Polynomial<type> getpoly(){return p;}
};