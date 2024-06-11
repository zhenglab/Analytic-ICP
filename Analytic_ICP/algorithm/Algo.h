#pragma once
/*!
*      \file Algo.h
*      \brief algorithm
*	   \author Wei Feng
*      \date 09/16/2020
*/
#ifndef MATHEMATICAL
#define MATHEMATICAL
#include<math.h>
#include <stdlib.h>
#include <time.h> 
#include "ICP\icpPointToPoint.h"

#define INFINITYNUM 16777216
#define MAX_VERTEX_NUM 32768

#define MAXITERCOUNT 10
#define PI 3.14159265358979323846264338327950288419716939937510
#define MAXNOMIALSIZE 64
#define MAXD 8
#define MAXDEGREE 32
#define MAXMATSIZE (MAXNOMIALSIZE*3)
typedef double coefmat[MAXMATSIZE];
typedef double double2[2];
typedef double double3[3];
typedef double double4[4];
typedef double double6[6];
typedef double double9[9];
typedef double double18[18];

#define RADTODEG(ang)			((ang) * 180.0 / PI)

inline double get_euclid_norm(VectorXd p)
{
	return sqrt(pow(p(0), 2) + pow(p(1), 2));
}

inline double get_euclid_dist(VectorXd p1, VectorXd p2)
{
	return sqrt(pow(p1(0) - p2(0), 2) + pow(p1(1) - p2(1), 2));
}

template<class T,int N>
T  __declspec(dllexport) get_euclid_dist(T* p1, T* p2)
{
	T dist = 0;
	for (int i = 0; i != N; ++i)
	{
		dist += pow(p1[i] - p2[i], 2);
	}
	dist = sqrt(dist);
	return dist;
}

inline void matrix2double3(double3 &d_set, Vector2d &p_set)
{
	d_set[0] = p_set(0);
	d_set[1] = p_set(1);
	d_set[2] = 1;
}

inline void matrix3double3(double3 &d_set, Vector3d &p_set)
{
	d_set[0] = p_set(0);
	d_set[1] = p_set(1);
	d_set[2] = p_set(2);
}

inline void matrix2double3(double3 *d_set, Vector2d *p_set, int num)
{
	for (int i = 0; i != num; ++i)
	{
		matrix2double3(d_set[i], p_set[i]);
	}
}

inline void matrix3double3(double3 *d_set, Vector3d *p_set, int num)
{
	for (int i = 0; i != num; ++i)
	{
		matrix3double3(d_set[i], p_set[i]);
	}
}


inline void double2matrix(Vector3d &p_set, double3 &d_set)
{
	for (int i = 0; i != 3; ++i)
	{
		p_set(i) = d_set[i];
	}
}

inline void double2matrix(Vector2d &p_set, double3 &d_set)
{
	for (int i = 0; i != 2; ++i)
	{
		p_set(i) = d_set[i];
	}
}


inline void double2matrix(Vector3d *p_set, double3 *d_set, int num)
{
	for (int i = 0; i != num; ++i)
	{
		double2matrix(p_set[i], d_set[i]);
	}
}

inline void double2matrix(Vector2d *p_set, double3 *d_set, int num)
{
	for (int i = 0; i != num; ++i)
	{
		double2matrix(p_set[i], d_set[i]);
	}
}


// (float,double) * (float,double) =  (double , float)
template<class TT1, class TT2, class TT3>
void  __declspec(dllexport) MultMatrix(TT1* M, TT2* M1, TT3* M2, int rows, int soms, int cols)
{
	//M - rows*soms    M1 - soms*cols   M2 - rows*cols
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j<cols; ++j)
		{
			M2[i*cols + j] = 0.0;
			for (int k = 0; k < soms; ++k)
			{
				M2[i*cols + j] += M[i*soms + k] * M1[k*cols + j];
			}
		}
	}
}

inline void GetNBB(MatrixXd &nbb, MatrixXd &matrixB)
{
	MatrixXd transposedB = matrixB.transpose();
	nbb = transposedB*matrixB;
}

inline void GetWu1(VectorXd &Wu1, MatrixXd &matrixB, VectorXd &l)
{
	MatrixXd transposedB = matrixB.transpose();
	Wu1 = transposedB*l;
}

//indirect adjustment model
inline void GetCorrection(VectorXd &correction, MatrixXd &matrixB, VectorXd &l, int N, int U)
{
	MatrixXd nbb(U, U);
	MatrixXd inverseForNbb(U, U);
	VectorXd Wu1(U);

	GetNBB(nbb, matrixB);
	inverseForNbb = nbb.inverse();
	GetWu1(Wu1, matrixB, l);
	correction = inverseForNbb*Wu1;
}



template<class T>
void __declspec(dllexport) AnalyticTrans3d_1(T(*mat_set)[MAXMATSIZE]
	, T(&srcObj)[3], T(&qstObj)[3])
{
	typedef T T3[3];
	qstObj[0] = mat_set[0][0];
	qstObj[1] = mat_set[0][1];
	qstObj[2] = mat_set[0][2];
	T3 tQst;

	T X[3] = { 0 };

	X[0] = srcObj[0];
	X[1] = srcObj[1];
	X[2] = srcObj[2];

	MultMatrix(mat_set[1], X, tQst, 3, 3, 1);
	qstObj[0] += tQst[0];
	qstObj[1] += tQst[1];
	qstObj[2] += tQst[2];
}

template<class T>
void __declspec(dllexport) Analyt3dTrans(T(*mat_set)[MAXMATSIZE]
	, int deg, T(&srcObj)[3], T(&qstObj)[3])
{
	typedef T T3[3];
	qstObj[0] = 0;
	qstObj[1] = 0;
	qstObj[2] = 0;
	T3 tQst;
	T coef[MAXNOMIALSIZE] = { 0 };
	int coef_size;
	for (int i = 0; i <= deg; ++i)
	{
		GetCoefOfTaylor(coef, coef_size, srcObj, i);

		MultMatrix(mat_set[i], coef, tQst, 3, coef_size, 1);
		qstObj[0] += tQst[0];
		qstObj[1] += tQst[1];
		qstObj[2] += tQst[2];
	}
}



template<class T>
void __declspec(dllexport) AnalyticTrans(T(*mat_set)[MAXMATSIZE]
	, T nomiVal[3], int deg
	, T(&srcObj)[3], T(&qstObj)[3])
{
	typedef T T2[2];
	qstObj[0] = 0;
	qstObj[1] = 0;
	qstObj[2] = 1;
	T2 tQst;
	for (int i = 0; i <= deg; ++i)
	{
		T *X = new T[i + 1];

		for (int j = 0; j <= i; ++j)
		{ 
			X[j] = Combination(i, j)
				/ Factorial(i)
				*pow(srcObj[0] - nomiVal[0], i - j)
				*pow(srcObj[1] - nomiVal[1], j);
		}

		MultMatrix(mat_set[i], X, tQst, 2, i + 1, 1);
		qstObj[0] += tQst[0];
		qstObj[1] += tQst[1];

		delete[]X;
	}
}


/*!
* the sum of n terms of the arithmetic sequence
*/
template<class T>
T __declspec(dllexport) Sn(T a0, T d, int n)
{
	T v = n*a0 + n*(n - 1)*d / 2;
	return v;
}

/*!
* the sum of n terms of the 1,3,6,10^ sequence
*/
inline int Sn136(int n)
{
	int v = n*(n + 1)*(n + 2) / 2;
	return v;
}


inline double Factorial(int n)
{
	if (n <= 1)
	{
		return 1.0;
	}
	else
	{
		return n*Factorial(n - 1);
	}
}

inline int Combination(int n, int k)
{
	if (n < k)
	{
		return 1;
	}

	int c = Factorial(k)*Factorial(n - k);
	return Factorial(n) / c;
}


/*!
* the point is moved by a randomly generated analytic map.
* \param v target
* \param mo the moved point
* \param num the point number
* \param deg the order of the Taylor series
*/
template<class TT1>
void __declspec(dllexport) AnalyticTransf(TT1(*v)[2], TT1(*mo)[2], int num, int deg)
{
	srand((unsigned)time(NULL));
	int sn = Sn(2, 2, deg + 1);
	double *pCoef = new double[sn];
	double nomiVal[2] = { 0 };

	FILE *fp = fopen("F:\\coef.csv", "w+");
	if (fp == NULL) {
		fprintf(stderr, "fopen() failed.\n");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i != sn; ++i)
	{
		pCoef[i] = 0.4 * rand() / RAND_MAX - 0.2;

		fprintf(fp, "%f\n", pCoef[i]);

	}

	nomiVal[0] = 0.1;
	nomiVal[1] = 0.1;

	fclose(fp);

	pCoef[2] = 1;
	pCoef[5] = 1;

	double *pTerm = new double[2];
	int sum = 0;

	for (int i = 0; i != num; ++i)
	{
		mo[i][0] = 0;
		mo[i][1] = 0;
		for (int j = 1; j <= deg; ++j)
		{
			sum = Sn(2, 2, j);
			for (int k = 0; k <= j; ++k)
			{
				pTerm[0] = pCoef[sum + k * 2 + 0]
					* Combination(j, k)
					*pow(v[i][0] - nomiVal[0], j - k)
					*pow(v[i][1] - nomiVal[1], k) / Factorial(j);

				pTerm[1] = pCoef[sum + k * 2 + 1]
					* Combination(j, k)
					*pow(v[i][0] - nomiVal[0], j - k)
					*pow(v[i][1] - nomiVal[1], k) / Factorial(j);

 				mo[i][0] += pTerm[0];
				mo[i][1] += pTerm[1];
			}
		}

	}

	delete[]pTerm;
	delete[]pCoef;
}


/*!
* the 3d point is moved by a randomly generated analytic map.
* \param v target
* \param mo the moved point
* \param num the point number
* \param deg the order of the Taylor series
*/

inline void AnalyticT3d(double(*v)[3], double(*mo)[3], int num, int deg)
{
	srand((unsigned)time(NULL));
	coefmat *mat_set = new coefmat[MAXDEGREE];
	memset(mat_set, 0, sizeof(coefmat)*MAXDEGREE);

	FILE *fp = fopen("F:\\coef.csv", "w+");
	if (fp == NULL) {
		fprintf(stderr, "fopen() failed.\n");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i <= deg; ++i)
	{
		for (int j = 0; j != MAXMATSIZE; ++j)
		{
			mat_set[i][j] = 0.2 * rand() / RAND_MAX - 0.1;
			fprintf(fp, "%f\n", mat_set[i][j]);
		}
	}

	fclose(fp);

	mat_set[1][0] = 1;
	mat_set[1][4] = 1;
	mat_set[1][8] = 1;
	double *pTerm = new double[3];
	int sum = 0;
	int item_num = 0;

	for (int i = 0; i != num; ++i)
	{
		mo[i][0] = mat_set[0][0];
		mo[i][1] = mat_set[0][1];
		mo[i][2] = mat_set[0][2];

		Analyt3dTrans(mat_set, deg, v[i], mo[i]);
	}
}

template<class T, int D>
void __declspec(dllexport) ExpectandVari(T ctr[D], T &vari, T(*ps)[D], int num)
{
	memset(ctr, 0, sizeof(T)*D);
	vari = 0;

	for (int i = 0; i != num; ++i)
	{
		for (int j = 0; j != D; ++j)
		{
			ctr[j] += ps[i][j] / num;
		}
	}

	for (int i = 0; i != num; ++i)
	{
		for (int j = 0; j != D; ++j)
		{
			T v = abs(ps[i][j] - ctr[j]);

			if (vari < v)
			{
				vari = v;
			}
		}
	}
}


inline bool GetCoefOfTaylor(double coef[MAXNOMIALSIZE], int &size, double p[3], int order)
{
	size = Combination(order + 2, 2);

	int index = 0;

	int i = order, j, k;
	double temp,temp1;
	while (1)
	{
		temp = 1;
		temp *= pow(p[0], i);
		temp /= Factorial(i);
		
		j = order - i;

		for (int k = 0; k <= j; ++k)
		{
			temp1 = temp;
			temp1 = temp1* pow(p[1], j - k);
			temp1 = temp1*pow(p[2], k);
			temp1 = temp1 / Factorial(k);
			temp1 = temp1 / Factorial(j - k);
			coef[index] = temp1;
			index++;
		}

		if (0 == i)
		{
			break;
		}

		i--;
	}

	if (index == size)
	{
		return true;
	}
	else
	{
		return false;
	}
}


inline void BFGS(MatrixXd &Uk, MatrixXd &Hk, VectorXd &delta_gk, VectorXd &delta_xk)
{
	MatrixXd coef = (delta_gk.transpose()*Hk*delta_gk) / (delta_xk.transpose()*delta_gk);

	MatrixXd div = delta_xk.transpose()*delta_gk;

	MatrixXd v = (delta_xk*delta_xk.transpose()) / div(0);

	MatrixXd l = (Hk*delta_gk*delta_xk.transpose()
		+ (Hk*delta_gk*delta_xk.transpose()).transpose()) / div(0);

	Uk = (1 + coef(0))*v - l;
}

#endif