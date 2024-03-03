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
#define MAXD 64
typedef double coefmat[MAXD];
typedef double double2[2];
typedef double double3[3];
typedef double double4[4];
typedef double double6[6];
typedef double double9[9];
typedef double double18[18];

#define RADTODEG(ang)			((ang) * 180.0 / PI)

typedef struct Point
{
	Point()
	{
		x = 0;
		y = 0;
	}
	int x;
	int y;
};

typedef struct Point2f
{
	Point2f()
	{
		x = 0;
		y = 0;
	}
	float x;
	float y;
};

typedef struct Rect
{
	Rect()
	{
		x = 0;
		y = 0;
		width = 0;
		height = 0;
	}
	int x;
	int y;
	int width;
	int height;
};


inline double get_euclid_dist(Point p, Rect r)
{
	return sqrt(pow(p.x - r.x, 2) + pow(p.y - r.y, 2));
}

inline double get_euclid_dist(CPoint2 p, Rect r)
{
	CPoint2 c;
	c[0] = r.x + r.width / 2;
	c[1] = r.y + r.height / 2;
	return sqrt(pow(p[0] - c[0], 2) + pow(p[1] - c[1], 2));
}

inline double get_euclid_norm(CPoint2 p)
{
	return sqrt(pow(p[0], 2) + pow(p[1], 2));
}

inline double get_euclid_dist(CPoint2 p1,CPoint2 p2)
{
	return sqrt(pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2));
}
inline double get_euclid_dist(CPoint2 p1, Point p2)
{
	return sqrt(pow(p1[0] - p2.x, 2) + pow(p1[1] - p2.y, 2));
}
inline double get_euclid_dist(Point2f p1, Point2f p2)
{
	return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

inline double get_euclid_dist(Point p1, Point p2)
{
	return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
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

inline void matrix2double3(double3 &d_set, CPoint2 &p_set)
{
	d_set[0] = p_set[0];
	d_set[1] = p_set[1];
	d_set[2] = 1;
}

inline void matrix3double3(double3 &d_set, CPoint3 &p_set)
{
	d_set[0] = p_set[0];
	d_set[1] = p_set[1];
	d_set[2] = p_set[2];
}

inline void matrix2double3(double3 *d_set, CPoint2 *p_set, int num)
{
	for (int i = 0; i != num; ++i)
	{
		matrix2double3(d_set[i], p_set[i]);
	}
}

inline void matrix3double3(double3 *d_set, CPoint3 *p_set, int num)
{
	for (int i = 0; i != num; ++i)
	{
		matrix3double3(d_set[i], p_set[i]);
	}
}

inline void double2matrix(CPoint2 &p_set, double3 &d_set)
{
	p_set[0] = d_set[0];
	p_set[1] = d_set[1];
}

inline void double2matrix(CPoint2 *p_set, double3 *d_set, int num)
{
	for (int i = 0; i != num; ++i)
	{
		double2matrix(p_set[i], d_set[i]);
	}
}


inline void double2matrix(CPoint3 &p_set, double3 &d_set)
{
	p_set[0] = d_set[0];
	p_set[1] = d_set[1];
	p_set[2] = d_set[2];
}

inline void double2matrix(CPoint3 *p_set, double3 *d_set, int num)
{
	for (int i = 0; i != num; ++i)
	{
		double2matrix(p_set[i], d_set[i]);
	}
}

inline double pointinregion(CPoint2 p, Rect r)
{
	CPoint2 r_center;
	r_center[0] = r.x + r.width / 2;
	r_center[1] = r.y + r.height / 2;
	return get_euclid_dist(p, r_center);
}

template<class TT1>
void  __declspec(dllexport) VectorMinus(TT1 *M, TT1 *M1, TT1 *M2)
{
	M2[0] = M[0] - M1[0];
	M2[1] = M[1] - M1[1];
	M2[2] = M[2] - M1[2];
}

//matrix  M2=M-M1   (double , float) - (double , float) = (double , float)
template<class TT1, class TT2, class TT3>
void  __declspec(dllexport) MatrixMinus(TT1* M, TT2* M1, TT3* M2, int nrows, int ncols)
{//M - nrows*ncols    M1 - nrows*ncols   M2 - nrows*ncols
	for (int i = 0; i < nrows; ++i)
	{
		for (int j = 0; j < ncols; ++j)
		{
			M2[i*ncols + j] = M[i*ncols + j] - M1[i*ncols + j];
		}
	}
}

// (float,double) * (float,double) =  (double , float)
template<class TT1, class TT2, class TT3>
void  __declspec(dllexport) MultMatrix(TT1* M, TT2* M1, TT3* M2, int rows, int cols)
{//M - rows*cols    M1 - cols*rows   M2 - rows*rows
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j<rows; ++j)
		{
			M2[i*rows + j] = 0.0;
			for (int k = 0; k < cols; ++k)
			{
				M2[i*rows + j] += M[i*cols + k] * M1[k*rows + j];
			}
		}
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


//matrix[3][4] change to  matrix0[4][3]
template<class T, class T0>
void __declspec(dllexport) MatrixTranspose(T* matrix, T0* matrix0, int nrows, int ncols)
{// matrix - nrows*ncols    matrix0 - ncols*nrows
	for (int i = 0; i < nrows; ++i)
	{
		for (int j = 0; j < ncols; ++j)
		{
			matrix0[j*nrows + i] = matrix[i*ncols + j];
		}
	}
}

/// <summary>   
/// Generating identity matrix
/// </summary>     
/// <param name="matrix"></param> 
/// <param name="nrows"></param>
template<class T>
void __declspec(dllexport) MatrixIdentity(T* matrix, int nrows)
{ //  matrix -  nrows*nrows
	memset(matrix, 0, nrows*nrows * sizeof(T));
	for (int i = 0; i < nrows; ++i)
	{
		matrix[i*nrows + i] = 1;
	}
}


/// <summary>   
/// QR factorization
/// </summary>  
template<class T>
int Bmaqr(T Myq[], T Myr[], T Mya[], int Mym, int Myn)
{
	int i, j, k, l, nn, p, jj;
	T u, alpha, w, t;

	if (Mym<Myn)
	{
		return(0);
	}

	for (i = 0; i <= Mym - 1; i++)
	{
		for (j = 0; j <= Mym - 1; j++)
		{
			l = i*Mym + j;
			Myq[l] = 0.0;
			if (i == j)
			{
				Myq[l] = 1.0;
			}

		}
	}

	for (i = 0; i <= Mym - 1; i++)
	{
		for (j = 0; j <= Myn - 1; j++)
		{
			l = i*Mym + j;
			Myr[l] = Mya[l];

		}
	}

	for (i = 0; i <= Mym - 2; i++)
	{
		for (j = i + 1; j <= Myn - 1; j++)
		{
			p = i*Mym + j;
			l = j*Mym + i;
			t = Myr[p];
			Myr[p] = Myr[l];
			Myr[l] = t;

		}
	}


	nn = Myn;
	if (Mym == Myn)
	{
		nn = Mym - 1;
	}
	for (k = 0; k <= nn - 1; k++)
	{
		u = 0.0;
		l = k*Myn + k;
		for (i = k; i <= Mym - 1; i++)
		{
			w = fabs(Myr[i*Myn + k]);
			if (w>u)
			{
				u = w;
			}

		}
		alpha = 0.0;
		for (i = k; i <= Mym - 1; i++)
		{
			t = Myr[i*Myn + k] / u;
			alpha = alpha + t*t;
		}
		if (Myr[l]>0.0)
		{
			u = -u;
		}
		alpha = u*sqrt(alpha);
		if (fabs(alpha) + 1.0 == 1.0)
		{
			return(0);
		}
		u = sqrt(2.0*alpha*(alpha - Myr[l]));
		if ((u + 1.0) != 1.0)
		{
			Myr[l] = (Myr[l] - alpha) / u;
			for (i = k + 1; i <= Mym - 1; i++)
			{
				p = i*Myn + k;
				Myr[p] = Myr[p] / u;
			}
			for (j = 0; j <= Mym - 1; j++)
			{
				t = 0.0;
				for (jj = k; jj <= Mym - 1; jj++)
				{
					t = t + Myr[jj*Myn + k] * Myq[jj*Mym + j];
				}
				for (i = k; i <= Mym - 1; i++)
				{
					p = i*Mym + j;
					Myq[p] = Myq[p] - 2.0*t*Myr[i*Myn + k];
				}

			}
			for (j = k + 1; j <= Myn - 1; j++)
			{
				t = 0.0;
				for (jj = k; jj <= Mym - 1; jj++)
				{
					t = t + Myr[jj*Myn + k] * Myr[jj*Myn + j];
				}
				for (i = k; i <= Mym - 1; i++)
				{
					p = i*Myn + j;
					Myr[p] = Myr[p] - 2.0*t*Myr[i*Myn + k];
				}

			}
			Myr[l] = alpha;
			for (i = k + 1; i <= Mym - 1; i++)
			{
				Myr[i*Myn + k] = 0.0;
			}

		}

	}

	for (i = 0; i <= Mym - 2; i++)
	{
		for (j = i + 1; j <= Myn - 1; j++)
		{
			p = i*Mym + j;
			l = j*Mym + i;
			t = Myr[p];
			Myr[p] = Myr[l];
			Myr[l] = t;

		}
	}
	return (1);

}


//Doolittle max elements solve equations 
//Doolittle Column principal elimination
template<class TT1, class TT2, class TT3>//>
void __declspec(dllexport) Doolittle(TT1* aa, TT2* bb, TT3* xx, int rows)
{// aa * xx = bb        root - xx[rows][rows]
	int k, i, j, t, ik;
	int* M = new int[rows];
	double  *s, *l, *u, *a, *b;
	double temp, smax = 0, *y, *x;
	s = new double[rows];
	l = new double[rows*rows];
	u = new double[rows*rows];
	a = new double[rows*rows];
	b = new double[rows];
	y = new double[rows];
	x = new double[rows];
	//  QA  =  LU
	for (i = 0; i<rows; ++i)
	{
		M[i] = 0;
		for (j = 0; j<rows; ++j)
		{
			a[i*rows + j] = aa[i*rows + j];
		}
	}
	for (k = 0; k<rows; ++k)
	{
		for (i = k; i<rows; ++i)
		{
			s[i] = a[i*rows + k];
			for (t = 0; t < k; ++t)
			{
				s[i] -= l[i*rows + t] * u[t*rows + k];
			}

			if (i == k)
			{
				smax = s[i];
				ik = i;
			}
			if (fabs(smax)<fabs(s[i]))
			{
				smax = s[i];
				ik = i;
			}
		}
		M[k] = ik;
		if (ik != k)
		{
			for (t = 0; t<k; ++t)
			{
				temp = l[k*rows + t];
				l[k*rows + t] = l[ik*rows + t];
				l[ik*rows + t] = temp;
			}
			for (t = k; t<rows; ++t)
			{
				temp = a[k*rows + t];
				a[k*rows + t] = a[ik*rows + t];
				a[ik*rows + t] = temp;
			}
			temp = s[k];
			s[k] = s[ik];
			s[ik] = temp;
		}
		u[k*rows + k] = s[k];
		if (k<rows - 1)
		{
			for (j = k + 1; j<rows; ++j)
			{
				u[k*rows + j] = a[k*rows + j];
				for (t = 0; t < k; ++t)
				{
					u[k*rows + j] -= l[k*rows + t] * u[t*rows + j];
				}
			}
			for (i = k + 1; i < rows; ++i)
			{
				l[i*rows + k] = s[i] / (u[k*rows + k] + 0.00001);
			}
		}
	}
	//Qb  =  Ly   AND   Ux  =   y
	for (j = 0; j<rows; ++j)
	{
		for (i = 0; i < rows; ++i)
		{
			b[i] = bb[i*rows + j];
		}
		for (k = 0; k<rows - 1; ++k)
		{
			t = M[k];
			temp = b[k];
			b[k] = b[t];
			b[t] = temp;
		}
		y[0] = b[0];
		for (i = 1; i<rows; ++i)
		{
			y[i] = b[i];
			for (t = 0; t<i; ++t)
			{
				y[i] -= l[i*rows + t] * y[t];
			}
		}
		x[rows - 1] = y[rows - 1] / (u[rows*rows - 1] + 0.00001);
		for (i = rows - 2; i>-1; --i)
		{
			x[i] = y[i];
			for (t = i + 1; t < rows; ++t)
			{
				x[i] -= u[i*rows + t] * x[t];
			}

			x[i] /= (u[i*rows + i] + 0.00001);
		}
		for (i = 0; i<rows; ++i)
		{
			xx[i*rows + j] = x[i];
		}
	}
	delete[]M;
	delete[]s;
	delete[]l;
	delete[]u;
	delete[]a;
	delete[]b;
	delete[]y;
	delete[]x;
	M = NULL;
	s = NULL;
	l = NULL;
	u = NULL;
	a = NULL;
	b = NULL;
	y = NULL;
	x = NULL;
}


/// <summary>   
/// matrix inversion
/// </summary>     
/// <param name="Matrix"></param> 
/// <param name="MatrixA">inverse matrix</param>
/// <param name="rows"></param>
//////////////////////////////////////////////////////////////////////////
// invertible matrice
template<class TT1, class TT2>
void __declspec(dllexport) MatrixAnti(TT1* Matrix, TT2* MatrixA, int rows)
{//  Matrix * MatrixA = I          I = E
	double* E = new double[rows*rows];
	MatrixIdentity(E, rows);

	//Doolittle solution
	Doolittle(Matrix, E, MatrixA, rows);
	delete[]E;
	E = NULL;
}

template<class TT1>
void __declspec(dllexport) GetNBB(TT1 nbb[], const TT1 *matrixB, int N, int U)
{
	TT1 *transposedB = new TT1[N*U];
	MatrixTranspose(matrixB, transposedB, N, U);
	MultMatrix(transposedB, matrixB, nbb, U, N, U);

	delete[] transposedB;
}

template<class TT1>
void __declspec(dllexport) GetNcc(TT1 Ncc[], const TT1 matrixA[], const TT1 nbb[], int S, int U)
{
	TT1 *inverseForNbb = new TT1[U*U];
	TT1 *transposedA = new TT1[S*U];
	TT1 *transit = new TT1[S*U];
	MatrixTranspose(matrixA, transposedA, S, U);
	MatrixAnti(nbb, inverseForNbb, U);
	MultMatrix(matrixA, inverseForNbb, transit, S, U, U);
	MultMatrix(transit, transposedA, Ncc, S, U, S);

	delete[] inverseForNbb;
	delete[] transposedA;
	delete[] transit;
}

template<class TT1>
void __declspec(dllexport) GetWu1(TT1 Wu1[], const TT1 *matrixB, const TT1 *l, int N, int U)
{
	TT1 *transposedB = new TT1[N*U];
	MatrixTranspose(matrixB, transposedB, N, U);
	MultMatrix(transposedB, l, Wu1, U, N, 1);

	delete[] transposedB;
}


//Indirect Adjustment with Constraints (in the case of singular matrices)
template<class TT1>
void __declspec(dllexport) GetCorrectionWithCondition4singular(TT1 correction[], const TT1 matrixA[], const TT1 w[], const TT1 *matrixB
	, const TT1 *l, int N, int U, int S)
{
	TT1 *nbb = new TT1[U*U];
	TT1 *Wu1 = new TT1[U];
	TT1 *transposedA = new TT1[S*U];
	TT1 *coefficient = new TT1[pow(U + S, 2)];
	TT1 *inverse4Coefficient = new TT1[pow(U + S, 2)];
	TT1 *Wxl = new TT1[U + S];
	TT1 *result = new TT1[U + S];
	GetNBB<TT1>(nbb, matrixB, N, U);
	GetWu1<TT1>(Wu1, matrixB, l, N, U);
	MatrixTranspose(matrixA, transposedA, S, U);

	memset(coefficient, 0, sizeof(TT1)*pow(U + S, 2));

	for (int i = 0; i != U; ++i)
	{
		TT1*p_nbb = nbb + i*U;
		TT1 *p_coe = coefficient + i*(U + S);
		for (int j = 0; j != U; ++j)
		{
			p_coe[j] = p_nbb[j];
		}
	}

	for (int i = 0; i != U; ++i)
	{
		TT1 *p_ta = transposedA + i*S;
		TT1 *p_coe = coefficient + i*(U + S);
		for (int j = 0; j != S; ++j)
		{
			p_coe[j + U] = p_ta[j];
		}
	}

	for (int i = 0; i != S; ++i)
	{
		TT1 *p_a = matrixA + i*U;
		TT1 *p_coe = coefficient + (i + U)*(U + S);
		for (int j = 0; j != U; ++j)
		{
			p_coe[j] = p_a[j];
		}
	}

	memcpy(Wxl, W, sizeof(TT1)*S);
	memcpy(Wxl + S, Wu1, sizeof(TT1)*U);

	MatrixAnti(coefficient, inverse4Coefficient, U + S);
	MultMatrix(inverse4Coefficient, Wxl, result, U + S, U + S);

	memcpy(correction, result, sizeof(TT1)*U);

	delete[] nbb;
	delete[] Wu1;
	delete[] transposedA;
	delete[] coefficient;
	delete[] inverse4Coefficient;
	delete[] Wxl;
	delete[] result;
}


//Indirect Adjustment with Constraints
template<class TT1>
void __declspec(dllexport) GetCorrectionWithCondition(TT1 correction[]
	, const TT1 matrixA[], const TT1 w[], const TT1 *matrixB
	, const TT1 *l, int N, int U, int S)
{
	TT1 *nbb = new TT1[U*U];
	TT1 *inverseForNbb = new TT1[U*U];
	TT1 *transposedA = new TT1[S*U];
	TT1 *Ncc = new TT1[S*S];
	TT1 *inverseForNcc = new TT1[S*S];
	TT1 *Wu1 = new TT1[U];

	GetNBB<TT1>(nbb, matrixB, N, U);
	MatrixAnti(nbb, inverseForNbb, U);
	MatrixTranspose(matrixA, transposedA, S, U);
	GetNcc<TT1>(Ncc, matrixA, nbb, S, U);
	MatrixAnti(Ncc, inverseForNcc, S);
	GetWu1<TT1>(Wu1, matrixB, l, N, U);

	TT1 *transit1 = new TT1[S*U];

	MultMatrix(inverseForNbb, transposedA, transit1, U, U, S);

	TT1 *transit2 = new TT1[S*U];

	MultMatrix(transit1, inverseForNcc, transit2, U, S, S);

	TT1 *transit3 = new TT1[U*U];

	MultMatrix(transit2, matrixA, transit3, U, S, U);

	TT1 *transit4 = new TT1[U*U];

	MultMatrix(transit3, inverseForNbb, transit4, U, U);

	TT1 *transit5 = new TT1[U*U];

	MatrixMinus(inverseForNbb, transit4, transit5, U, U);

	TT1 *transit6 = new TT1[U];

	MultMatrix(transit5, Wu1, transit6, U, U, 1);

	TT1 *transit7 = new TT1[S*U];

	MultMatrix(inverseForNbb, transposedA, transit7, U, U, S);

	TT1 *transit8 = new TT1[S*U];

	MultMatrix(transit7, inverseForNcc, transit8, U, S, S);

	TT1 *transit9 = new TT1[U];

	MultMatrix(transit8, w, transit9, U, S, 1);

	MatrixMinus(transit6, transit9, correction, U, 1);

	delete[] nbb;
	delete[] inverseForNbb;
	delete[] transposedA;
	delete[] Ncc;
	delete[] inverseForNcc;
	delete[] Wu1;
	delete[] transit1;
	delete[] transit2;
	delete[] transit3;
	delete[] transit4;
	delete[] transit5;
	delete[] transit6;
	delete[] transit7;
	delete[] transit8;
	delete[] transit9;
}

/// <summary>   
/// Conditional adjustment with restricted conditions
/// </summary>     
/// <param name="correction"></param> 
/// <param name="MatrixA"></param>
/// <param name="MatrixB"></param>
/// <param name="MatrixC"></param>
/// <param name="W"></param>
/// <param name="Wx"></param>
/// <param name="C"></param>
/// <param name="N"></param>
/// <param name="U"></param>
/// <param name="S"></param>
//////////////////////////////////////////////////////////////////////////
template<class T>
void __declspec(dllexport) GetConditionCorrectionWithCondition(T correction[], const T matrixA[], const T matrixB[], const T matrixC[], const T W[], const T Wx[], int C, int N, int U, int S)
{
	T *transposedA = new T[N*C];
	T *transposedB = new T[U*C];
	T *transposedC = new T[U*S];
	T *Naa = new T[C*C];
	T *inverseForNaa = new T[C*C];
	T *Nbb = new T[U*U];
	T *inverseForNbb = new T[U*U];
	T *Ncc = new T[S*S];
	T *inverseForNcc = new T[S*S];
	T *We = new T[U];
	T *x = new T[U];
	MatrixTranspose(matrixA, transposedA, C, N);
	MultMatrix(matrixA, transposedA, Naa, C, N, C);
	MatrixAnti(Naa, inverseForNaa, C);
	MatrixTranspose(matrixB, transposedB, C, U);
	T *transit1 = new T[U*C];
	MultMatrix(transposedB, inverseForNaa, transit1, U, C, C);
	MultMatrix(transit1, matrixB, Nbb, U, C, U);
	MatrixAnti(Nbb, inverseForNbb, U);

	T *transit2 = new T[U*C];
	MultMatrix(transposedB, inverseForNaa, transit2, U, C, C);
	MultMatrix(transit2, W, We, U, C, 1);

	T *transit3 = new T[S*U];
	MultMatrix(matrixC, inverseForNbb, transit3, S, U, U);
	MatrixTranspose(matrixC, transposedC, S, U);
	MultMatrix(transit3, transposedC, Ncc, S, U, S);
	MatrixAnti(Ncc, inverseForNcc, S);


	T *transit4 = new T[U*S];
	MultMatrix(inverseForNbb, transposedC, transit4, U, U, S);
	T *transit5 = new T[U*S];
	MultMatrix(transit4, inverseForNcc, transit5, U, S, S);

	T *transit6 = new T[U*U];
	MultMatrix(transit5, matrixC, transit6, U, S, U);
	T *transit7 = new T[U*U];
	MultMatrix(transit6, inverseForNbb, transit7, U, U);
	MatrixMinus(inverseForNbb, transit7, transit7, U, U);
	T *transit8 = new T[U];
	MultMatrix(transit7, We, transit8, U, U, 1);

	T *transit9 = new T[U];
	MultMatrix(transit5, Wx, transit9, U, S, 1);

	MatrixMinus(transit9, transit9, correction, U, 1);
	MatrixMinus(correction, transit8, correction, U, 1);
	MatrixMinus(correction, transit9, correction, U, 1);


	delete[] transposedA;
	delete[] transposedB;
	delete[] transposedC;
	delete[] Naa;
	delete[] inverseForNaa;
	delete[] Nbb;
	delete[] inverseForNbb;
	delete[] Ncc;
	delete[] inverseForNcc;
	delete[] We;
	delete[] x;
	delete[] transit1;
	delete[] transit2;
	delete[] transit3;
	delete[] transit4;
	delete[] transit5;
	delete[] transit6;
	delete[] transit7;
	delete[] transit8;
	delete[] transit9;
}


//indirect adjustment model
template<class TT1>
void __declspec(dllexport) GetCorrection(TT1 correction[], const TT1 *matrixB, const TT1 *l, int N, int U)
{
	TT1 *nbb = new TT1[U*U];
	TT1 *inverseForNbb = new TT1[U*U];
	TT1 *Wu1 = new TT1[U];

	GetNBB<TT1>(nbb, matrixB, N, U);
	MatrixAnti(nbb, inverseForNbb, U);
	GetWu1<TT1>(Wu1, matrixB, l, N, U);

	MultMatrix(inverseForNbb, Wu1, correction, U, U, 1);

	delete[] nbb;
	delete[] inverseForNbb;
	delete[] Wu1;
}

template<class TT1>
TT1 __declspec(dllexport) GetVecNorm(TT1* point0, int D)
{
	TT1 result = 0;

	for (int i = 0; i<D; ++i)
	{
		result += pow(point0[i], 2);
	}

	return sqrt(result);
}

template<class T>
void __declspec(dllexport) Power2Trans(T(&linearT)[9], T(&o2)[9], T(&srcObj)[3], T(&qstObj)[3])
{
	typedef T T3[3];
	T3 tObj,o2term;

	MultMatrix(linearT, srcObj, qstObj, 3, 3, 1);

	o2term[0] = pow(srcObj[0], 2);
	o2term[1] = srcObj[0] * srcObj[1];
	o2term[2] = pow(srcObj[1], 2);

	MultMatrix(o2, o2term, tObj, 2, 3, 1);

	qstObj[0] += tObj[0];
	qstObj[1] += tObj[1];
}

template<class T>
void __declspec(dllexport) AnalyticTrans3d_1(T(*mat_set)[MAXD]
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
void __declspec(dllexport) AnalyticTrans3d_2(T(&mat_set)[MAXD]
	, T(&srcObj)[3], T(&qstObj)[3])
{
	typedef T T3[3];
	T X[1024] = { 0 };
	X[0] = srcObj[0] * srcObj[0] / 2;
	X[1] = srcObj[0] * srcObj[1];
	X[2] = srcObj[0] * srcObj[2];
	X[3] = srcObj[1] * srcObj[1] / 2;
	X[4] = srcObj[1] * srcObj[2];
	X[5] = srcObj[2] * srcObj[2] / 2;

	MultMatrix(mat_set, X, qstObj, 3, 6, 1);
}



template<class T>
void __declspec(dllexport) AnalyticTrans(T(*mat_set)[MAXD], int deg
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
				/Factorial(i)
				*pow(srcObj[0], i - j)
				*pow(srcObj[1], j);
		}

		MultMatrix(mat_set[i], X, tQst, 2, i + 1, 1);
		qstObj[0] += tQst[0];
		qstObj[1] += tQst[1];

		delete[]X;
	}
}

template<class T>
void __declspec(dllexport) Power3Trans(T(&linearT)[9], T(&o2)[9], T(&o3)[9]
	, T(&srcObj)[3], T(&qstObj)[3])
{
	typedef T T2[2];
	typedef T T3[3];
	typedef T T4[4];
	T2 tObj;
	T3 o2term;
	T4 o3term;

	MultMatrix(linearT, srcObj, qstObj, 3, 3, 1);

	o2term[0] = pow(srcObj[0], 2);
	o2term[1] = srcObj[0] * srcObj[1];
	o2term[2] = pow(srcObj[1], 2);

	MultMatrix(o2, o2term, tObj, 2, 3, 1);

	qstObj[0] += tObj[0];
	qstObj[1] += tObj[1];

	o3term[0] = pow(srcObj[0], 3);
	o3term[1] = pow(srcObj[0], 2) * srcObj[1];
	o3term[2] = srcObj[0] * pow(srcObj[1], 2);
	o3term[3] = pow(srcObj[1], 3);

	MultMatrix(o3, o3term, tObj, 2, 4, 1);

	qstObj[0] += tObj[0];
	qstObj[1] += tObj[1];
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
					*pow(v[i][0], j - k)
					*pow(v[i][1], k) / Factorial(j);

				pTerm[1] = pCoef[sum + k * 2 + 1]
					* Combination(j, k)
					*pow(v[i][0], j - k)
					*pow(v[i][1], k) / Factorial(j);

 				mo[i][0] += pTerm[0];
				mo[i][1] += pTerm[1];
			}
		}

	}
	
	delete[]pTerm;
}


/*!
* the 3d point is moved by a randomly rank 2 generated analytic map.
* \param v target
* \param mo the moved point
* \param num the point number
*/

inline void AnalyticT3d(double(*v)[3], double(*mo)[3], int num)
{
	srand((unsigned)time(NULL));
	int sn = Sn136(3);
	double pCoef[1024] = { 0 };

	FILE *fp = fopen("F:\\coef.csv", "w+");
	if (fp == NULL) {
		fprintf(stderr, "fopen() failed.\n");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i != sn; ++i)
	{
		pCoef[i] = 0.1 * rand() / RAND_MAX - 0.05;

		fprintf(fp, "%f\n", pCoef[i]);

	}

	fclose(fp);

	pCoef[3] = 1;
	pCoef[7] = 1;
	pCoef[11] = 1;
	double *pTerm = new double[3];
	int sum = 0;
	int item_num = 0;

	for (int i = 0; i != num; ++i)
	{
		mo[i][0] = pCoef[0];
		mo[i][1] = pCoef[1];
		mo[i][2] = pCoef[2];

		for (int j = 1; j <= 2; ++j)
		{
			sum = Sn136(j - 1) * 3 + 3;
			item_num = Combination(3 + j - 1, 3 - 1);
			double X[1024] = {0};

			if (1 == j)
			{
				X[0] = v[i][0];
				X[1] = v[i][1];
				X[2] = v[i][2];
			}
			else if (2 == j)
			{
				X[0] = v[i][0] * v[i][0] / 2;
				X[1] = v[i][0] * v[i][1];
				X[2] = v[i][0] * v[i][2];
				X[3] = v[i][1] * v[i][1] / 2;
				X[4] = v[i][1] * v[i][2];
				X[5] = v[i][2] * v[i][2] / 2;
			}
			//int x_num = 0;
			//for (int k = 0; k <= j; ++k)
			//{
			//	X[x_num] = pow(v[i][0], j - k) / Factorial(j - k);
			//	for (int l = 0; l <= k; ++l)
			//	{
			//		X[x_num] *= pow(v[i][1], k - l) / Factorial(k - l);
			//		X[x_num] *= pow(v[i][2], l) / Factorial(l);
			//		x_num++;
			//	}
			//}

			MultMatrix(pCoef + sum, X, pTerm, 3, item_num, 1);
			mo[i][0] += pTerm[0];
			mo[i][1] += pTerm[1];
			mo[i][2] += pTerm[2];
		}

	}
}

#endif