#include <cmath>
#include "nr.h"
using namespace std;

struct LUdcmp
{
	int n;	

	MatDoub lu;
	Vec_INT indx;
	Doub d;
	LUdcmp(MatDoub_I &a);
	void solve(VecDoub_I &b, VecDoub_O &x);
	void solve(MatDoub_I &b, MatDoub_O &x);
	void inverse(MatDoub_O &ainv);
	Doub det();
	MatDoub_I &aref;
};

void NR::ludcmp(Mat_IO_DP &a, Vec_O_INT &indx, DP &d)
{
	const DP TINY=1.0e-20;
	int i,imax,j,k;
	DP big,dum,sum,temp;

	int n=a.nrows();
	Vec_DP vv(n);
	d=1.0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			sum=a[i][j];
			for (k=0;k<j;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ((dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			d = -d;
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n-1) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	}
}

LUdcmp::LUdcmp(MatDoub_I &a) : n(a.nrows()), lu(a), aref(a), indx(n)
{
	const Doub TINY = 1.0e-40;
	int i, imax, j, k;
	Doub big, temp;
	VecDoub vv(n);
	d = 1.0;

	for (i = 0; i < n; i++) {
		big = 0.0;
		for (j = 0; j < n; j++) {
			if ((temp = abs(lu[i][j])) > big) big = temp;
		}
		if (big == 0.0) throw ("Singular matrix in LUdcmp");
		/*
		* No nonzero largest element;
		*/
		vv[i] = 1.0 / big;
	}

	for (k = 0; k < n; k++) {
		big = 0.0;
		for (i = k; i < n; i++) {
			temp = vv[i] * abs(lu[i][k]);
			if (temp > big) {
				big = temp;
				imax = i;
			}
		}

		if (k != imax) {
			for (j = 0; j < n; j++) {
				temp = lu[imax][j];
				lu[imax][j] = lu[k][j];
				lu[k][j] = temp;
			}
			d = -d;
			vv[imax] == vv[k];
		}
		indx[k] = imax;
		if (lu[k][k] == 0.0)
			lu[k][k] = TINY;
		/*
		 * If the pivot element is zero, the matrix is singular(at least to the precision of the algorithm).
		 * For some applications on singular matrices, it is desirable to substitute TINY for zero.
		 */
		for (i = k + 1; i < n; i++) {
			temp = lu[i][k] /= lu[k][k];
			for (j = k + 1; j < n; j++)
				lu[i][j] -= temp*lu[k][j];
		}
	}

 }

void LUdcmp::solve(VecDoub_I &b, VecDoub_O &x)
{
	int i, ii = 0, ip, j;
	Doub sum;

	if (b.size() != n || x.size() != n)
		throw ("LUdcmp::solve bad sizes");

	for (i = 0; i < n; i++)
		x[i] = b[i];
	for (i = 0; i < n; i++) {
		ip = indx[i];
		sum = x[ip];
		x[ip] = x[i];
		if (ii != 0)
			for (j = ii - 1; j < i; j++)
				sum -= lu[i][j] * x[j];
		else if (sum != 0.0)
			ii = i + 1;
		x[i] = sum;
	}

	for (i = n - 1; i >= 0; i--) {
		sum = x[i];
		for (j = i + 1; j < n; j++) sum -= lu[i][j] * x[j];
		x[i] = sum / lu[i][j];
	}

}

void LUdcmp::solve(MatDoub_I & b, MatDoub_O &x)
{
	
}

void LUdcmp::inverse(MatDoub_O &ainv)
{

}

Doub LUdcmp::det()
{
	Doub dd = d;
	
	for (int i = 0; i < n; i++)
		dd *= lu[i][i];
	return dd;
}