#include <cmath>
#include "nr.h"
#include "iomanip"
using namespace std;

void NR::gaussj(Mat_IO_DP &a, Mat_IO_DP &b)
{
	int i,icol,irow,j,k,l,ll, z, ls;
	DP big,dum,pivinv;

	int n=a.nrows();
	int m=b.ncols();
	
	Vec_INT indxc(n),indxr(n),ipiv(n);

	for (j=0;j<n;j++) ipiv[j]=0;	//These integer arrays are used for bookkeeping on the pivoting
	for (i=0;i<n;i++) {				//This is the main loop over the columns to be reduced.
		big=0.0;
		cout << endl << "i =" << i << endl;
		for (j=0;j<n;j++)			//This is the outer loop of the search for a pivot element;
			if (ipiv[j] != 1)
				for (k=0;k<n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
							cout << "irow = " << irow << setw(12) << "icol = " << icol << endl;
						}
					}
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l]);
			for (l=0;l<m;l++) SWAP(b[irow][l],b[icol][l]);

			cout << endl << "change Matrix a : " << endl;
			for (z = 0; z<n; z++) {
				for (ls = 0; ls<n; ls++) cout << setw(12) << a[z][ls];
				cout << endl;
			}

			cout << endl << "change vector b : " << endl;
			for (z = 0; z<n; z++) {
				for (ls = 0; ls<m; ls++) cout << setw(12) << b[z][ls];
				cout << endl;
			}	

		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix");
		pivinv=1.0/a[icol][icol];
		cout << "pivinv = " << pivinv << endl;
		a[icol][icol]=1.0;
		for (l=0;l<n;l++) a[icol][l] *= pivinv;
		for (l=0;l<m;l++) b[icol][l] *= pivinv;

		cout << endl << "pivinv Matrix a : " << endl;
		for (z = 0; z<n; z++) {
			for (ls = 0; ls<n; ls++) cout << setw(12) << a[z][ls];
			cout << endl;
		}
		cout << endl << "pivinv vector b : " << endl;
		for (z = 0; z<n; z++) {
			for (ls = 0; ls<m; ls++) cout << setw(12) << b[z][ls];
			cout << endl;
		}
		cout << "Current icol = " << icol << setw(20) <<  "Currenct irow = " << irow << endl;
		for (ll=0;ll<n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				cout << "dum = " << dum << endl;
				a[ll][icol]=0.0;
				cout << "a[" << ll << "][" << icol << "] = " << a[ll][icol] << endl;
				for (l = 0; l < n; l++) {
					a[ll][l] -= a[icol][l] * dum;
					cout << "a[" << ll << "][" << l << "] -= a[" << icol << "][" << l << "] * " << dum << ";" << endl;
				}
				cout << "a[" << ll << "][" << icol << "] = " << a[ll][icol] << endl;
				for (l = 0; l < m; l++) {
					b[ll][l] -= b[icol][l] * dum;
					cout << "b[" << ll << "][" << l << "] -= b[" << icol << "][" << l << "] * " << dum << ";" << endl;
				}
			}

		/* Out put the change*/
		cout << endl << "Matrix a : " << i << endl;
		for (z = 0; z<n; z++) {
			for (ls = 0; ls<n; ls++) cout << setw(12) << a[z][ls];
			cout << endl;
		}
		cout << endl << "vector b : " << i <<  endl;
		for (z = 0; z<n; z++) {
			for (ls = 0; ls<m; ls++) cout << setw(12) << b[z][ls];
			cout << endl;
		}
	}
	for (l=n-1;l>=0;l--) {
		if (indxr[l] != indxc[l])
			for (k=0;k<n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
}
