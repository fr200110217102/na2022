#include<bits/stdc++.h>
using namespace std;

int main(){
    double UFL = 1.0 / 2;
    double OFL = 2 * (2 - 1.0 / 4);
    cout << "UFL = " << UFL << endl;
    cout << "OFL = " << OFL << endl;
    vector <double> num;
    vector <double> subnum;
    num.push_back(0);
    double base = 1.0 / 2;
    for (int e = -1; e <= 1; ++ e) {
        double M;
        for (int i = 0; i < 2; ++ i)
            for (int j = 0; j < 2; ++ j) {
                for (int k = 0; k < 2; ++ k) {
                    M = i + j / 2.0 + k / 4.0;
                    if (i > 0) {
                        num.push_back(M * base);
                        num.push_back(-M * base);
                    }
                    else {
                        if (e == -1) {
                            subnum.push_back(M * base);
                            subnum.push_back(-M * base);
                        }
                    }
                }
            }
        base *= 2;
    }
    sort(num.begin(), num.end());
    cout << "numbers(" << num.size() << "):\n";
    for (double x : num) cout << x << ",\n"[x == num.back()];
    sort(subnum.begin(), subnum.end());
    cout << "subnumbers(" << subnum.size() << "):\n";
    for (double x : subnum) cout << x << ",\n"[x == subnum.back()];
}