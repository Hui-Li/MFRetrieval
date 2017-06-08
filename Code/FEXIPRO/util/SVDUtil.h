#ifndef SVDUTIL_H
#define SVDUTIL_H

#include "Base.h"
#include <mlpack/core.hpp>

using namespace arma;

int SVD(const mat &P_t, const int m, const int n, Matrix &u, Matrix &v, const double SIGMA){

	mat U_t;
	vec s;
	mat V;

	// see: http://arma.sourceforge.net/docs.html#svd_econ
//	svd_econ(U_t, s, V, P_t, "both", "std");
	svd_econ(U_t, s, V, P_t);

	U_t = U_t.t();

	double *uData = new double[m * m];
	double *vData = new double[m * n];

	for (int rowIndex = 0; rowIndex < m; rowIndex++) {
		for (int colIndex = 0; colIndex < m; colIndex++) {
			uData[rowIndex * m + colIndex] = s[rowIndex] * U_t(rowIndex, colIndex);
		}

	}

	vector<double> sum(m);
	sum[0] = s[0];
	for (int colIndex = 1; colIndex < m; colIndex++) {
		sum[colIndex] = sum[colIndex - 1] + s[colIndex];
	}

	int checkDim = 0;
	for (int colIndex = 0; colIndex < m; colIndex++) {
		if(sum[colIndex] / sum[m - 1] >= SIGMA) {
			checkDim = colIndex;
			break;
		}
	}

	for(int rowIndex = 0; rowIndex < n; rowIndex++){
		for (int colIndex = 0; colIndex < m; colIndex++) {
			vData[rowIndex * m + colIndex] = V(rowIndex, colIndex);
		}
	}

	u.init(uData, m, m);
	v.init(vData, n, m);

	return checkDim;
}

int SVD(const mat &P_t, const int m, const int n, Matrix &u, Matrix &v, vector<double> &sigmaValues, const double SIGMA){

	mat U_t;
	vec s;
	mat V;

	// see: http://arma.sourceforge.net/docs.html#svd_econ
//	svd_econ(U_t, s, V, P_t, "both", "std");
	svd_econ(U_t, s, V, P_t);

	U_t = U_t.t();

	double *uData = new double[m * m];
	double *vData = new double[m * n];

	for (int rowIndex = 0; rowIndex < m; rowIndex++) {
		for (int colIndex = 0; colIndex < m; colIndex++) {
			uData[rowIndex * m + colIndex] = s[rowIndex] * U_t(rowIndex, colIndex);
		}

	}

	sigmaValues.resize(m);
	vector<double> sum(m);
	sigmaValues[0] = s[0] / s[m-1];
	sum[0] = s[0];
	for (int colIndex = 1; colIndex < m; colIndex++) {
		sigmaValues[colIndex] = s[colIndex] / s[m-1];
		sum[colIndex] = sum[colIndex - 1] + s[colIndex];
	}

	int checkDim = 0;
	for (int colIndex = 0; colIndex < m; colIndex++) {
		if(sum[colIndex] / sum[m - 1] >= SIGMA) {
			checkDim = colIndex;
			break;
		}
	}

	for(int rowIndex = 0; rowIndex < n; rowIndex++){
		for (int colIndex = 0; colIndex < m; colIndex++) {
			vData[rowIndex * m + colIndex] = V(rowIndex, colIndex);
		}
	}

	u.init(uData, m, m);
	v.init(vData, n, m);

	return checkDim;
}

void SVD(const mat &P_t, const int m, const int n, Matrix &u, Matrix &v, double *vSubNorms) {

	mat U_t;
	vec s;
	mat V;

	// see: http://arma.sourceforge.net/docs.html#svd_econ
//	svd_econ(U_t, s, V, P_t, "both", "std");
	svd_econ(U_t, s, V, P_t);

	U_t = U_t.t();

	u.init(m, m);
	v.init(n, m);

	for (int rowIndex = 0; rowIndex < m; rowIndex++) {
		double *uPtr = u.getRowPtr(rowIndex);
		for (int colIndex = 0; colIndex < m; colIndex++) {
			uPtr[colIndex] = s[rowIndex] * U_t(rowIndex, colIndex);
		}
	}

	for(int rowIndex = 0; rowIndex < n; rowIndex++){
		double *subVNorm = &vSubNorms[rowIndex * m];
		double norm = 0;
		double *vPtr = v.getRowPtr(rowIndex);
		for (int colIndex = m-1; colIndex >= 0; colIndex--) {
			vPtr[colIndex] = V(rowIndex, colIndex);
			norm += V(rowIndex, colIndex) * V(rowIndex, colIndex);
			subVNorm[colIndex] = sqrt(norm);
		}
	}

}
#endif //SVDUTIL_H
