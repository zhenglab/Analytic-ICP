#ifndef FITTING
#define FITTING
/*!
*      \file Fitting.h
*      \brief Fitting
*      \author Wei Feng
*      \date 10/08/2020
*
*/
#include "Algo.h"

#define MAXNEWICPITERCOUNT 64

template<class T>
void  __declspec(dllexport) GetMatB4Plane(T matB[], T(*point_set)[3], int point_num)
{
	memset(matB, 0, sizeof(T) * 3 * point_num);
	for (int i = 0; i != point_num; ++i)
	{
		matB[3 * i + 0] = point_set[i][0];
		matB[3 * i + 1] = point_set[i][1];
		matB[3 * i + 2] = point_set[i][2];
	}
}

template<class T>
void  __declspec(dllexport) GetMatB4Affine(T matB[], T(*point0)[3], T(*point)[3], int point_num)
{
	memset(matB, 0, sizeof(T) * 2 * 6 * point_num);

	for (int i = 0; i != point_num; ++i)
	{
		matB[i * 2 * 6 + 6 * 0 + 0] = point0[i][0];
		matB[i * 2 * 6 + 6 * 0 + 1] = point0[i][1];
		matB[i * 2 * 6 + 6 * 0 + 2] = 0;
		matB[i * 2 * 6 + 6 * 0 + 3] = 0;
		matB[i * 2 * 6 + 6 * 0 + 4] = 1;
		matB[i * 2 * 6 + 6 * 0 + 5] = 0;

		matB[i * 2 * 6 + 6 * 1 + 0] = 0;
		matB[i * 2 * 6 + 6 * 1 + 1] = 0;
		matB[i * 2 * 6 + 6 * 1 + 2] = point0[i][0];
		matB[i * 2 * 6 + 6 * 1 + 3] = point0[i][1];
		matB[i * 2 * 6 + 6 * 1 + 4] = 0;
		matB[i * 2 * 6 + 6 * 1 + 5] = 1;
	}
}

template<class T>
void  __declspec(dllexport) GetMatB3d4AffineQR(T matB[], T(&ortMat)[9]
	, T(*point0)[3], int point_num)
{
	memset(matB, 0, sizeof(T) * 3 * 9 * point_num);

	for (int i = 0; i != point_num; ++i)
	{
		matB[i * 3 * 9 + 9 * 0 + 0] = 1;
		matB[i * 3 * 9 + 9 * 0 + 1] = 0;
		matB[i * 3 * 9 + 9 * 0 + 2] = 0;
		matB[i * 3 * 9 + 9 * 0 + 3] = ortMat[0] * point0[i][0];
		matB[i * 3 * 9 + 9 * 0 + 4] = ortMat[0] * point0[i][1];
		matB[i * 3 * 9 + 9 * 0 + 5] = ortMat[0] * point0[i][2];
		matB[i * 3 * 9 + 9 * 0 + 6] = ortMat[1] * point0[i][1];
		matB[i * 3 * 9 + 9 * 0 + 7] = ortMat[1] * point0[i][2];
		matB[i * 3 * 9 + 9 * 0 + 8] = ortMat[2] * point0[i][2];

		matB[i * 3 * 9 + 9 * 1 + 0] = 0;
		matB[i * 3 * 9 + 9 * 1 + 1] = 1;
		matB[i * 3 * 9 + 9 * 1 + 2] = 0;
		matB[i * 3 * 9 + 9 * 1 + 3] = ortMat[3] * point0[i][0];
		matB[i * 3 * 9 + 9 * 1 + 4] = ortMat[3] * point0[i][1];
		matB[i * 3 * 9 + 9 * 1 + 5] = ortMat[3] * point0[i][2];
		matB[i * 3 * 9 + 9 * 1 + 6] = ortMat[4] * point0[i][1];
		matB[i * 3 * 9 + 9 * 1 + 7] = ortMat[4] * point0[i][2];
		matB[i * 3 * 9 + 9 * 1 + 8] = ortMat[5] * point0[i][2];

		matB[i * 3 * 9 + 9 * 2 + 0] = 0;
		matB[i * 3 * 9 + 9 * 2 + 1] = 0;
		matB[i * 3 * 9 + 9 * 2 + 2] = 1;
		matB[i * 3 * 9 + 9 * 2 + 3] = ortMat[6] * point0[i][0];
		matB[i * 3 * 9 + 9 * 2 + 4] = ortMat[6] * point0[i][1];
		matB[i * 3 * 9 + 9 * 2 + 5] = ortMat[6] * point0[i][2];
		matB[i * 3 * 9 + 9 * 2 + 6] = ortMat[7] * point0[i][1];
		matB[i * 3 * 9 + 9 * 2 + 7] = ortMat[7] * point0[i][2];
		matB[i * 3 * 9 + 9 * 2 + 8] = ortMat[8] * point0[i][2];
	}
}


template<class T>
void  __declspec(dllexport) GetMatB4AffineQR(T matB[], T(&ortMat)[9], T(*point0)[3], int point_num)
{
	memset(matB, 0, sizeof(T) * 2 * 5 * point_num);

	for (int i = 0; i != point_num; ++i)
	{
		matB[i * 2 * 5 + 5 * 0 + 0] = 1;
		matB[i * 2 * 5 + 5 * 0 + 1] = 0;
		matB[i * 2 * 5 + 5 * 0 + 2] = ortMat[0] * point0[i][0];
		matB[i * 2 * 5 + 5 * 0 + 3] = ortMat[0] * point0[i][1];
		matB[i * 2 * 5 + 5 * 0 + 4] = ortMat[1] * point0[i][1];

		matB[i * 2 * 5 + 5 * 1 + 0] = 0;
		matB[i * 2 * 5 + 5 * 1 + 1] = 1;
		matB[i * 2 * 5 + 5 * 1 + 2] = ortMat[3] * point0[i][0];
		matB[i * 2 * 5 + 5 * 1 + 3] = ortMat[3] * point0[i][1];
		matB[i * 2 * 5 + 5 * 1 + 4] = ortMat[4] * point0[i][1];
	}
}

template<class T>
void  __declspec(dllexport) GetMatB4ADMM(T matB[], int deg
	, T(*point0)[3], int point_num)
{
	int i, j;
	double comb_num;
	int c = 2 * (deg + 1);
	double one_fac = 1.0 / Factorial(deg);
	for (i = 0; i != point_num; ++i)
	{
		for (j = 0; j <= deg; ++j)
		{
			comb_num = Combination(deg, j);
			matB[i * 2 * c + c * 0 + j] = one_fac*comb_num
				*pow(point0[i][0], deg - j)*pow(point0[i][1], j);
			matB[i * 2 * c + c * 1 + j] = 0;

			matB[i * 2 * c + c * 0 + deg + 1 + j] = 0;
			matB[i * 2 * c + c * 1 + deg + 1 + j] = one_fac*comb_num
				*pow(point0[i][0], deg - j)*pow(point0[i][1], j);
		}
	}
}


template<class T>
void  __declspec(dllexport) GetMatB4Analytic(T matB[], int deg, int c
	, T(*point0)[3], int point_num)
{
	int e_i;
	int i, j, k;
	double one_fac;
	int comb_num;
	for (i = 0; i != point_num; ++i)
	{
		e_i = 0;
		for (j = 0; j <= deg; ++j)
		{
			one_fac = 1.0 / Factorial(j);
			for (k = 0; k <= j; ++k)
			{
				comb_num = Combination(j, k);
				matB[i * 2 * c + c * 0 + e_i + k] = one_fac*comb_num
					*pow(point0[i][0], j - k)*pow(point0[i][1], k);
				matB[i * 2 * c + c * 1 + e_i + k] = 0;

				matB[i * 2 * c + c * 0 + e_i + j + 1 + k] = 0;
				matB[i * 2 * c + c * 1 + e_i + j + 1 + k] = one_fac*comb_num
					*pow(point0[i][0], j - k)*pow(point0[i][1], k);
			}
			e_i += 2 * (j + 1);
		}
	}
}

template<class T>
void  __declspec(dllexport) GetMatB4Cubic(T matB[], T(*point0)[3], int point_num)
{
	memset(matB, 0, sizeof(T) * 2 * 20 * point_num);
	for (int i = 0; i != point_num; ++i)
	{
		matB[i * 2 * 20 + 20 * 0 + 0] = 1;
		matB[i * 2 * 20 + 20 * 0 + 1] = 0;
		matB[i * 2 * 20 + 20 * 0 + 2] = point0[i][0];
		matB[i * 2 * 20 + 20 * 0 + 3] = point0[i][1];
		matB[i * 2 * 20 + 20 * 0 + 4] = 0;
		matB[i * 2 * 20 + 20 * 0 + 5] = 0;
		matB[i * 2 * 20 + 20 * 0 + 6] = pow(point0[i][0], 2);
		matB[i * 2 * 20 + 20 * 0 + 7] = point0[i][0] * point0[i][1];
		matB[i * 2 * 20 + 20 * 0 + 8] = pow(point0[i][1], 2);
		matB[i * 2 * 20 + 20 * 0 + 9] = 0;
		matB[i * 2 * 20 + 20 * 0 + 10] = 0;
		matB[i * 2 * 20 + 20 * 0 + 11] = 0;
		matB[i * 2 * 20 + 20 * 0 + 12] = pow(point0[i][0], 3);
		matB[i * 2 * 20 + 20 * 0 + 13] = pow(point0[i][0], 2) * point0[i][1];
		matB[i * 2 * 20 + 20 * 0 + 14] = point0[i][0] * pow(point0[i][1], 2);
		matB[i * 2 * 20 + 20 * 0 + 15] = pow(point0[i][1], 3);
		matB[i * 2 * 20 + 20 * 0 + 16] = 0;
		matB[i * 2 * 20 + 20 * 0 + 17] = 0;
		matB[i * 2 * 20 + 20 * 0 + 18] = 0;
		matB[i * 2 * 20 + 20 * 0 + 19] = 0;

		matB[i * 2 * 20 + 20 * 1 + 0] = 0;
		matB[i * 2 * 20 + 20 * 1 + 1] = 1;
		matB[i * 2 * 20 + 20 * 1 + 2] = 0;
		matB[i * 2 * 20 + 20 * 1 + 3] = 0;
		matB[i * 2 * 20 + 20 * 1 + 4] = point0[i][0];
		matB[i * 2 * 20 + 20 * 1 + 5] = point0[i][1];
		matB[i * 2 * 20 + 20 * 1 + 6] = 0;
		matB[i * 2 * 20 + 20 * 1 + 7] = 0;
		matB[i * 2 * 20 + 20 * 1 + 8] = 0;
		matB[i * 2 * 20 + 20 * 1 + 9] = pow(point0[i][0], 2);
		matB[i * 2 * 20 + 20 * 1 + 10] = point0[i][0] * point0[i][1];
		matB[i * 2 * 20 + 20 * 1 + 11] = pow(point0[i][1], 2);
		matB[i * 2 * 20 + 20 * 1 + 12] = 0;
		matB[i * 2 * 20 + 20 * 1 + 13] = 0;
		matB[i * 2 * 20 + 20 * 1 + 14] = 0;
		matB[i * 2 * 20 + 20 * 1 + 15] = 0;
		matB[i * 2 * 20 + 20 * 1 + 16] = pow(point0[i][0], 3);
		matB[i * 2 * 20 + 20 * 1 + 17] = pow(point0[i][0], 2) * point0[i][1];
		matB[i * 2 * 20 + 20 * 1 + 18] = point0[i][0] * pow(point0[i][1], 2);
		matB[i * 2 * 20 + 20 * 1 + 19] = pow(point0[i][1], 3);
	}
}


template<class T>
void  __declspec(dllexport) GetMatB4Quad3d(T matB[], T(*point0)[3], int point_num)
{
	memset(matB, 0, sizeof(T) * 3 * 18 * point_num);

	for (int i = 0; i != point_num; ++i)
	{
		matB[i * 3 * 18 + 18 * 0 + 0] = 0.5*pow(point0[i][0], 2);
		matB[i * 3 * 18 + 18 * 0 + 1] = 0;
		matB[i * 3 * 18 + 18 * 0 + 2] = 0;
		matB[i * 3 * 18 + 18 * 0 + 3] = point0[i][0] * point0[i][1];
		matB[i * 3 * 18 + 18 * 0 + 4] = 0;
		matB[i * 3 * 18 + 18 * 0 + 5] = 0;
		matB[i * 3 * 18 + 18 * 0 + 6] = point0[i][0] * point0[i][2];
		matB[i * 3 * 18 + 18 * 0 + 7] = 0;

		matB[i * 3 * 18 + 18 * 0 + 8] = 0;
		matB[i * 3 * 18 + 18 * 0 + 9] = 0.5*pow(point0[i][1], 2);
		matB[i * 3 * 18 + 18 * 0 + 10] = 0;
		matB[i * 3 * 18 + 18 * 0 + 11] = 0;
		matB[i * 3 * 18 + 18 * 0 + 12] = point0[i][1] * point0[i][2];
		matB[i * 3 * 18 + 18 * 0 + 13] = 0;
		matB[i * 3 * 18 + 18 * 0 + 14] = 0;
		matB[i * 3 * 18 + 18 * 0 + 15] = 0.5*pow(point0[i][2], 2);
		matB[i * 3 * 18 + 18 * 0 + 16] = 0;
		matB[i * 3 * 18 + 18 * 0 + 17] = 0;


		matB[i * 3 * 18 + 18 * 1 + 0] = 0;
		matB[i * 3 * 18 + 18 * 1 + 1] = 0.5*pow(point0[i][0], 2);
		matB[i * 3 * 18 + 18 * 1 + 2] = 0;
		matB[i * 3 * 18 + 18 * 1 + 3] = 0;
		matB[i * 3 * 18 + 18 * 1 + 4] = point0[i][0] * point0[i][1];
		matB[i * 3 * 18 + 18 * 1 + 5] = 0;
		matB[i * 3 * 18 + 18 * 1 + 6] = 0;
		matB[i * 3 * 18 + 18 * 1 + 7] = point0[i][0] * point0[i][2];

		matB[i * 3 * 18 + 18 * 1 + 8] = 0;
		matB[i * 3 * 18 + 18 * 1 + 9] = 0;
		matB[i * 3 * 18 + 18 * 1 + 10] = 0.5*pow(point0[i][1], 2);
		matB[i * 3 * 18 + 18 * 1 + 11] = 0;
		matB[i * 3 * 18 + 18 * 1 + 12] = 0;
		matB[i * 3 * 18 + 18 * 1 + 13] = point0[i][1] * point0[i][2];
		matB[i * 3 * 18 + 18 * 1 + 14] = 0;
		matB[i * 3 * 18 + 18 * 1 + 15] = 0;
		matB[i * 3 * 18 + 18 * 1 + 16] = 0.5*pow(point0[i][2], 2);
		matB[i * 3 * 18 + 18 * 1 + 17] = 0;


		matB[i * 3 * 18 + 18 * 2 + 0] = 0;
		matB[i * 3 * 18 + 18 * 2 + 1] = 0;
		matB[i * 3 * 18 + 18 * 2 + 2] = 0.5*pow(point0[i][0], 2);
		matB[i * 3 * 18 + 18 * 2 + 3] = 0;
		matB[i * 3 * 18 + 18 * 2 + 4] = 0;
		matB[i * 3 * 18 + 18 * 2 + 5] = point0[i][0] * point0[i][1];
		matB[i * 3 * 18 + 18 * 2 + 6] = 0;
		matB[i * 3 * 18 + 18 * 2 + 7] = 0;

		matB[i * 3 * 18 + 18 * 2 + 8] = point0[i][0] * point0[i][2];
		matB[i * 3 * 18 + 18 * 2 + 9] = 0;
		matB[i * 3 * 18 + 18 * 2 + 10] = 0;
		matB[i * 3 * 18 + 18 * 2 + 11] = 0.5*pow(point0[i][1], 2);
		matB[i * 3 * 18 + 18 * 2 + 12] = 0;
		matB[i * 3 * 18 + 18 * 2 + 13] = 0;
		matB[i * 3 * 18 + 18 * 2 + 14] = point0[i][1] * point0[i][2];
		matB[i * 3 * 18 + 18 * 2 + 15] = 0;
		matB[i * 3 * 18 + 18 * 2 + 16] = 0;
		matB[i * 3 * 18 + 18 * 2 + 17] = 0.5*pow(point0[i][2], 2);
	}
}





template<class T>
void  __declspec(dllexport) GetMatB4Quadratic(T matB[], T(*point0)[3], int point_num)
{
	memset(matB, 0, sizeof(T) * 2 * 12 * point_num);
	for (int i = 0; i != point_num; ++i)
	{
		matB[i * 2 * 12 + 12 * 0 + 0] = 1;
		matB[i * 2 * 12 + 12 * 0 + 1] = 0;
		matB[i * 2 * 12 + 12 * 0 + 2] = point0[i][0];
		matB[i * 2 * 12 + 12 * 0 + 3] = point0[i][1];
		matB[i * 2 * 12 + 12 * 0 + 4] = 0;
		matB[i * 2 * 12 + 12 * 0 + 5] = 0;
		matB[i * 2 * 12 + 12 * 0 + 6] = pow(point0[i][0], 2);
		matB[i * 2 * 12 + 12 * 0 + 7] = point0[i][0]*point0[i][1];
		matB[i * 2 * 12 + 12 * 0 + 8] = pow(point0[i][1], 2);;
		matB[i * 2 * 12 + 12 * 0 + 9] = 0;
		matB[i * 2 * 12 + 12 * 0 + 10] = 0;
		matB[i * 2 * 12 + 12 * 0 + 11] = 0;

		matB[i * 2 * 12 + 12 * 1 + 0] = 0;
		matB[i * 2 * 12 + 12 * 1 + 1] = 1;
		matB[i * 2 * 12 + 12 * 1 + 2] = 0;
		matB[i * 2 * 12 + 12 * 1 + 3] = 0;
		matB[i * 2 * 12 + 12 * 1 + 4] = point0[i][0];
		matB[i * 2 * 12 + 12 * 1 + 5] = point0[i][1];
		matB[i * 2 * 12 + 12 * 1 + 6] = 0;
		matB[i * 2 * 12 + 12 * 1 + 7] = 0;
		matB[i * 2 * 12 + 12 * 1 + 8] = 0;
		matB[i * 2 * 12 + 12 * 1 + 9] = pow(point0[i][0], 2);
		matB[i * 2 * 12 + 12 * 1 + 10] = point0[i][0] * point0[i][1];
		matB[i * 2 * 12 + 12 * 1 + 11] = pow(point0[i][1], 2);
	}
}

template<class T>
void  __declspec(dllexport) GetMatA4Perspect(T matA[], T matP[], int point_num)
{
	memset(matA, 0, sizeof(T)*(8 + point_num));
	matA[0] = matP[4];
	matA[1] = -1 * matP[3];
	matA[3] = -1 * matP[1];
	matA[4] = matP[0];
}

template<class T>
void  __declspec(dllexport) GetMatB4Perspect(T matB[], T matP[], T(*point0)[3], T(*point)[3], int point_num)
{
	memset(matB, 0, sizeof(T) * 2 * 2 * point_num);

	for (int i = 0; i != point_num; ++i)
	{
		matB[i * 2 * 2 + 2 * 0 + 0] = point[i][0] * point0[i][0];
		matB[i * 2 * 2 + 2 * 0 + 1] = point[i][0] * point0[i][1];


		matB[i * 2 * 2 + 2 * 1 + 0] = point[i][1] * point0[i][0];
		matB[i * 2 * 2 + 2 * 1 + 1] = point[i][1] * point0[i][1];
	}
}

template<class T>
void  __declspec(dllexport) GetL4Plane(T l[], T param[], T(*point_set)[3], int point_num)
{
	memset(l, 0, sizeof(T)*point_num);

	for (int i = 0; i != point_num; ++i)
	{
		l[i] = 1 - param[0] * point_set[i][0] - param[1] * point_set[i][1] - param[2] * point_set[i][2];
	}
}

template<class T>
void  __declspec(dllexport) GetL4Affine(T l[], T param[], T(*point0)[3], T(*point)[3], int point_num)
{
	memset(l, 0, sizeof(T)*point_num * 2);

	for (int i = 0; i != point_num; ++i)
	{
		l[2 * i + 0] = point[i][0] - param[0] * point0[i][0] - param[1] * point0[i][1] - param[4];
		l[2 * i + 1] = point[i][1] - param[2] * point0[i][0] - param[3] * point0[i][1] - param[5];
	}
}

template<class T>
void  __declspec(dllexport) GetL3d4AffineQR(T l[], T(&ortMat)[9], T(&tVec)[3]
	, T param[], T(*point0)[3], T(*point)[3], int point_num)
{
	memset(l, 0, sizeof(T)*point_num * 3);

	for (int i = 0; i != point_num; ++i)
	{
		l[3 * i + 0] = point[i][0] - tVec[0] - param[0]
			- ortMat[0] * param[3] * point0[i][0]
			- ortMat[0] * param[4] * point0[i][1]
			- ortMat[1] * param[6] * point0[i][1]
			- ortMat[0] * param[5] * point0[i][2]
			- ortMat[1] * param[7] * point0[i][2] 
			- ortMat[2] * param[8] * point0[i][2];
		l[3 * i + 1] = point[i][1] - tVec[1] - param[1]
			- ortMat[3] * param[3] * point0[i][0]
			- ortMat[3] * param[4] * point0[i][1]
			- ortMat[4] * param[6] * point0[i][1]
			- ortMat[3] * param[5] * point0[i][2]
			- ortMat[4] * param[7] * point0[i][2]
			- ortMat[5] * param[8] * point0[i][2];
		l[3 * i + 2] = point[i][2] - tVec[2] - param[2]
			- ortMat[6] * param[3] * point0[i][0]
			- ortMat[6] * param[4] * point0[i][1]
			- ortMat[7] * param[6] * point0[i][1]
			- ortMat[6] * param[5] * point0[i][2]
			- ortMat[7] * param[7] * point0[i][2]
			- ortMat[8] * param[8] * point0[i][2];
	}
}


template<class T>
void  __declspec(dllexport) GetL4AffineQR(T l[], T(&ortMat)[9]
	, T param[], T(*point0)[3], T(*point)[3], int point_num)
{
	memset(l, 0, sizeof(T)*point_num * 2);

	for (int i = 0; i != point_num; ++i)
	{
		l[2 * i + 0] = point[i][0] - ortMat[2] - param[0] 
			- ortMat[0] * param[2] * point0[i][0] 
			- ortMat[0] * param[3] * point0[i][1] 
			- ortMat[1] * param[4] * point0[i][1];
		l[2 * i + 1] = point[i][1] - ortMat[5] - param[1] 
			- ortMat[3] * param[2] * point0[i][0] 
			- ortMat[3] * param[3] * point0[i][1] 
			- ortMat[4] * param[4] * point0[i][1];
	}
}

template<class T>
void  __declspec(dllexport) GetL4Analytic(T l[], T(*mat_set)[MAXMATSIZE], int deg
	, T(*point0)[3], T(*point)[3], int point_num)
{
	memset(l, 0, sizeof(T)*point_num * 2);
	int i, j, k;
	double one_fac;
	int comb_num;
	for (i = 0; i != point_num; ++i)
	{
		l[2 * i + 0] = point[i][0];
		l[2 * i + 1] = point[i][1];

		for (j = 0; j <= deg; ++j)
		{
			one_fac = 1.0 / Factorial(j);
			for (k = 0; k <= j; ++k)
			{
				comb_num = Combination(j, k);
				l[2 * i + 0] -= one_fac*comb_num
					*pow(point0[i][0], j - k)*pow(point0[i][1], k)*mat_set[j][k];
				l[2 * i + 1] -= one_fac*comb_num
					*pow(point0[i][0], j - k)*pow(point0[i][1], k)*mat_set[j][j + k + 1];
			}
		}
	}
}
 
template<class T>
void  __declspec(dllexport) GetL4Cubic(T l[], T param[], T(*point0)[3], T(*point)[3], int point_num)
{
	memset(l, 0, sizeof(T)*point_num * 2);

	for (int i = 0; i != point_num; ++i)
	{
		l[2 * i + 0] = point[i][0] - param[0] - param[2] * point0[i][0] - param[3] * point0[i][1]
			- param[6] * pow(point0[i][0], 2) - param[7] * point0[i][0] * point0[i][1]
			- param[8] * pow(point0[i][1], 2) - param[12] * pow(point0[i][0], 3) 
			- param[13] * pow(point0[i][0], 2)* point0[i][1] 
			- param[14] * point0[i][0] * pow(point0[i][1], 2) 
			- param[15] * pow(point0[i][1], 3);

		l[2 * i + 1] = point[i][1] - param[1] - param[4] * point0[i][0] - param[5] * point0[i][1]
			- param[9] * pow(point0[i][0], 2) - param[10] * point0[i][0] * point0[i][1]
			- param[11] * pow(point0[i][1], 2) - param[16] * pow(point0[i][0], 3)
			- param[17] * pow(point0[i][0], 2)* point0[i][1]
			- param[18] * point0[i][0] * pow(point0[i][1], 2)
			- param[19] * pow(point0[i][1], 3);
	}
}


template<class T>
void  __declspec(dllexport) GetL4Quad3d(T l[]
	, T param[], T(*point0)[3], T(*point)[3], int point_num)
{
	memset(l, 0, sizeof(T)*point_num * 3);

	for (int i = 0; i != point_num; ++i)
	{
		l[3 * i + 0] = point[i][0] - 
			- 0.5*param[0] * pow(point0[i][0], 2)
			- param[3] * point0[i][0] * point0[i][1]
			- param[6] * point0[i][0] * point0[i][2]
			- 0.5*param[9] * pow(point0[i][1], 2)
			- param[12] * point0[i][1] * point0[i][2]
			- 0.5*param[15] * pow(point0[i][2], 2);

		l[3 * i + 1] = point[i][1] - 
			- 0.5*param[1] * pow(point0[i][0], 2)
			- param[4] * point0[i][0] * point0[i][1]
			- param[7] * point0[i][0] * point0[i][2]
			- 0.5*param[10] * pow(point0[i][1], 2)
			- param[13] * point0[i][1] * point0[i][2]
			- 0.5*param[16] * pow(point0[i][2], 2);

		l[3 * i + 2] = point[i][2] - 
			- 0.5*param[2] * pow(point0[i][0], 2)
			- param[5] * point0[i][0] * point0[i][1]
			- param[8] * point0[i][0] * point0[i][2]
			- 0.5*param[11] * pow(point0[i][1], 2)
			- param[14] * point0[i][1] * point0[i][2]
			- 0.5*param[17] * pow(point0[i][2], 2);
	}
}


template<class T>
void  __declspec(dllexport) GetL4Quadratic(T l[], T param[], T(*point0)[3], T(*point)[3], int point_num)
{
	memset(l, 0, sizeof(T)*point_num * 2);

	for (int i = 0; i != point_num; ++i)
	{
		l[2 * i + 0] = point[i][0] - param[0] - param[2] * point0[i][0] - param[3] * point0[i][1]
			- param[6] * pow(point0[i][0], 2) - param[7] * point0[i][0] * point0[i][1]
			- param[8] * pow(point0[i][1], 2);

		l[2 * i + 1] = point[i][1] - param[1] - param[4] * point0[i][0] - param[5] * point0[i][1]
			- param[9] * pow(point0[i][0], 2) - param[10] * point0[i][0] * point0[i][1]
			- param[11] * pow(point0[i][1], 2);
	}
}

template<class T>
void  __declspec(dllexport) GetL4Perspect(T l[], T matP[], T(*point0)[3], T(*point)[3], int point_num)
{
	memset(l, 0, sizeof(T)*point_num * 2);

	for (int i = 0; i != point_num; ++i)
	{
		l[2 * i + 0] = point0[i][0] * matP[0] + point0[i][1] * matP[1] + matP[2]
			- point[i][0] * point0[i][0] * matP[6]
			- point[i][0] * point0[i][1] * matP[7]
			- point[i][0];
		l[2 * i + 1] = point0[i][0] * matP[3] + point0[i][1] * matP[4] + matP[5]
			- point[i][1] * point0[i][0] * matP[6]
			- point[i][1] * point0[i][1] * matP[7]
			- point[i][1];
	}
}

template<class T>
void  __declspec(dllexport) GetDeltaVec4Plane(T DeltaV[], T param[]
	, T(*point_set)[3], int point_num)
{
	T *matB = new T[3 * point_num];
	T *l = new T[point_num];

	GetMatB4Plane(matB, point_set, point_num);
	GetL4Plane(l, param, point_set, point_num);
	GetCorrection(DeltaV, matB, l, point_num, 3);

	delete[] matB;
	delete[] l;
}


template<class T>
void  __declspec(dllexport) GetDeltaVec4Affine(T DeltaV[], T param[]
	, T(*point0)[3], T(*point)[3], int point_num)
{
	T *matB = new T[2 * point_num * 6];
	T *l = new T[point_num * 2];

	GetMatB4Affine(matB, point0, point, point_num);
	GetL4Affine(l, param, point0, point, point_num);
	GetCorrection(DeltaV, matB, l, point_num * 2, 6);

	delete[] matB;
	delete[] l;
}


template<class T>
void  __declspec(dllexport) GetDVec3d4AffineQR(T DVec[], T(&ortMat)[9], T(&tVec)[3]
	, T param[], T(*point0)[3], T(*point)[3], int point_num)
{
	T *matB = new T[3 * point_num * 9];
	T *l = new T[point_num * 3];

	GetMatB3d4AffineQR(matB, ortMat, point0, point_num);
	GetL3d4AffineQR(l, ortMat, tVec, param, point0, point, point_num);
	GetCorrection(DVec, matB, l, point_num * 3, 9);

	delete[] matB;
	delete[] l;
}

template<class T>
void  __declspec(dllexport) GetDeltaVec4AffineQR(T DeltaV[], T(&ortMat)[9], T param[]
	, T(*point0)[3], T(*point)[3], int point_num)
{
	T *matB = new T[2 * point_num * 5];
	T *l = new T[point_num * 2];

	GetMatB4AffineQR(matB, ortMat, point0, point_num);
	GetL4AffineQR(l, ortMat, param, point0, point, point_num);
	GetCorrection(DeltaV, matB, l, point_num * 2, 5);

	delete[] matB;
	delete[] l;
}

template<class T>
T  __declspec(dllexport) GetDeltaVec4ADMM(T(&delta)[MAXMATSIZE], T(*mat_set)[MAXMATSIZE]
	, int deg, T(*point0)[3], T(*point)[3], int point_num)
{
	int sn = 2 * (deg + 1);
	T *matB = new T[2 * point_num * sn];
	T *l = new T[point_num * 2];
	memset(matB, 0, sizeof(T) * 2 * sn*point_num);

	GetMatB4ADMM(matB, deg, point0, point_num);
	GetL4Analytic(l, mat_set, deg, point0, point, point_num);
	GetCorrection(delta, matB, l, point_num * 2, sn);

	T error = GetVecNorm(delta, sn);

	delete[] matB;
	delete[] l;

	return error;
}


template<class T>
T  __declspec(dllexport) GetDeltaVec4Analytic(T(*delta_set)[MAXMATSIZE], T(*mat_set)[MAXMATSIZE]
	, int deg, T(*point0)[3], T(*point)[3], int point_num)
{
	int sn = Sn(2, 2, deg + 1);
	T *matB = new T[2 * point_num * sn];
	T *l = new T[point_num * 2];
	T *deltaV = new T[sn];
	memset(matB, 0, sizeof(T) * 2 * sn*point_num);

	GetMatB4Analytic(matB, deg, sn, point0, point_num);
	GetL4Analytic(l, mat_set, deg, point0, point, point_num);
	GetCorrection(deltaV, matB, l, point_num * 2, sn);

	int sn_temp;
	for (int i = 0; i <= deg; ++i)
	{
		sn_temp = Sn(2, 2, i);
		memcpy(delta_set[i], deltaV + sn_temp, sizeof(T) * 2 * (i + 1));
	}

	T error = GetVecNorm(deltaV, sn);

	delete[] matB;
	delete[] l;
	delete[] deltaV;

	return error;
}

template<class T>
void  __declspec(dllexport) GetDeltaVec4Cubic(T DeltaV[], T param[]
	, T(*point0)[3], T(*pointAff)[3], int point_num)
{
	T *matB = new T[2 * point_num * 20];
	T *l = new T[point_num * 2];

	GetMatB4Cubic(matB, point0, point_num);
	GetL4Cubic(l, param, point0, pointAff, point_num);
	GetCorrection(DeltaV, matB, l, point_num * 2, 20);

	delete[] matB;
	delete[] l;
}

template<class T>
void  __declspec(dllexport) GetDeltaVec4Quad3d(T DeltaV[], T param[]
	, T(*point0)[3], T(*pointAff)[3], int point_num)
{
	T *matB = new T[3 * point_num * 18];
	T *l = new T[point_num * 3];

	GetMatB4Quad3d(matB, point0, point_num);
	GetL4Quad3d(l, param, point0, pointAff, point_num);
	GetCorrection(DeltaV, matB, l, point_num * 3, 18);

	delete[] matB;
	delete[] l;
}

template<class T>
void  __declspec(dllexport) GetDeltaVec4Quadratic(T DeltaV[], T param[]
	, T(*point0)[3], T(*pointAff)[3], int point_num)
{
	T *matB = new T[2 * point_num * 12];
	T *l = new T[point_num * 2];

	GetMatB4Quadratic(matB, point0, point_num);
	GetL4Quadratic(l, param, point0, pointAff, point_num);
	GetCorrection(DeltaV, matB, l, point_num * 2, 12);

	delete[] matB;
	delete[] l;
}

template<class T>
void  __declspec(dllexport) GetDeltaVec4Perspect(T DeltaV[], T matP[],
	T(*point0)[3], T(*point)[3], int point_num)
{
	T *matB = new T[2 * 2 * point_num];
	T *l = new T[point_num * 2];

	GetMatB4Perspect(matB, matP, point0, point, point_num);
	GetL4Perspect(l, matP, point0, point, point_num);

	GetCorrection(DeltaV, matB, l, point_num * 2, 2);

	delete[] matB;
	delete[] l;
}


template<class T>
void  __declspec(dllexport) PlaneFitting(T param[3], T(*point_set)[3], int point_num, double eps, int jt)
{
	T *deltaV = new T[3];

	param[0] = 0;
	param[1] = 0;
	param[2] = 1;

	for (int i = 0; i != jt; ++i)
	{
		GetDeltaVec4Plane(deltaV, param, point_set, point_num);

		param[0] += deltaV[0];
		param[1] += deltaV[1];
		param[2] += deltaV[2];

		if (GetVecNorm(deltaV, 3) < eps)
		{
			break;
		}
	}

	delete[] deltaV;
}


template<class T>
void  __declspec(dllexport) AffineFitting(T param[6], T(*point0)[3], T(*point)[3]
	, int point_num, double eps, int jt)
{
	T *deltaV = new T[6];

	param[0] = 1;
	param[1] = 0;
	param[2] = 0;
	param[3] = 1;
	param[4] = 0;
	param[5] = 0;

	for (int i = 0; i != jt; ++i)
	{
		GetDeltaVec4Affine(deltaV, param, point0, point, point_num);

		param[0] += deltaV[0];
		param[1] += deltaV[1];
		param[2] += deltaV[2];
		param[3] += deltaV[3];
		param[4] += deltaV[4];
		param[5] += deltaV[5];

		if (GetVecNorm(deltaV, 6) < eps)
		{
			break;
		}
	}

	delete[] deltaV;
}


template<class T>
void  __declspec(dllexport) Affine3dFittingQR(T param[9], T(&ortMat)[9], T(&tVec)[3]
	, T(*point0)[3], T(*point)[3], int point_num, double eps, int jt)
{
	T *deltaV = new T[9];

	param[0] = 0;
	param[1] = 0;
	param[2] = 0;
	param[3] = 1;
	param[4] = 0;
	param[5] = 0;
	param[6] = 1;
	param[7] = 0;
	param[8] = 1;

	for (int i = 0; i != jt; ++i)
	{
		GetDVec3d4AffineQR(deltaV, ortMat, tVec, param, point0, point, point_num);

		param[0] += deltaV[0];
		param[1] += deltaV[1];
		param[2] += deltaV[2];
		param[3] += deltaV[3];
		param[4] += deltaV[4];
		param[5] += deltaV[5];
		param[6] += deltaV[6];
		param[7] += deltaV[7];
		param[8] += deltaV[8];

		if (GetVecNorm(deltaV, 9) < eps)
		{
			break;
		}
	}

	delete[] deltaV;
}



template<class T>
void  __declspec(dllexport) AffineFittingQR(T param[5], T(&ortMat)[9]
	, T(*point0)[3], T(*point)[3], int point_num, double eps, int jt)
{
	T *deltaV = new T[5];

	param[0] = 0;
	param[1] = 0;
	param[2] = 1;
	param[3] = 0;
	param[4] = 1;

	for (int i = 0; i != jt; ++i)
	{
		GetDeltaVec4AffineQR(deltaV, ortMat, param, point0, point, point_num);

		param[0] += deltaV[0];
		param[1] += deltaV[1];
		param[2] += deltaV[2];
		param[3] += deltaV[3];
		param[4] += deltaV[4];

		if (GetVecNorm(deltaV, 5) < eps)
		{
			break;
		}
	}

	delete[] deltaV;
}

template<class T>
void  __declspec(dllexport) ADMMFitting(T(*mat_set)[MAXMATSIZE], int deg
	, T(*point0)[3], T(*point)[3], int point_num, double eps, int jt)
{
	typedef T Td[MAXMATSIZE];
	Td *delta_set = new Td[MAXDEGREE];
	memset(delta_set, 0, sizeof(Td)*MAXDEGREE);
	T error,tempE;
	for (int i = 0; i != jt; ++i)
	{
		error = 0;
		for (int j = 0; j <= deg; ++j)
		{
			tempE = GetDeltaVec4ADMM(delta_set[j], mat_set, j, point0, point, point_num);
			for (int k = 0; k < 2 * (j + 1); ++k)
			{
				mat_set[j][k] += delta_set[j][k];
			}
			error += tempE;
		}


		if (error < eps)
		{
			break;
		}
	}
}

template<class T>
void  __declspec(dllexport) AnalyticFitting(T(*mat_set)[MAXMATSIZE], int deg
	, T(*point0)[3], T(*pointAff)[3], int point_num, double eps, int jt)
{
	typedef T Td[MAXMATSIZE];
	Td *delta_set = new Td[MAXDEGREE];
	memset(delta_set, 0, sizeof(Td)*MAXDEGREE);
	for (int i = 0; i != jt; ++i)
	{
		T error = GetDeltaVec4Analytic(delta_set, mat_set, deg, point0, pointAff, point_num);

		for (int j = 0; j <= deg; ++j)
		{
			for (int k = 0; k < 2 * (j + 1); ++k)
			{
				mat_set[j][k] += delta_set[j][k];
			}
		}

		if (error < eps)
		{
			break;
		}
	}
}

template<class T>
void  __declspec(dllexport) CubicFitting(T param[20]
	, T(*point0)[3], T(*pointAff)[3], int point_num, double eps, int jt)
{
	T *deltaV = new T[20];
	for (int i = 0; i != jt; ++i)
	{
		getDeltaVec4Cubic(deltaV, param, point0, pointAff, point_num);

		param[0] += deltaV[0];
		param[1] += deltaV[1];
		param[2] += deltaV[2];
		param[3] += deltaV[3];
		param[4] += deltaV[4];
		param[5] += deltaV[5];
		param[6] += deltaV[6];
		param[7] += deltaV[7];
		param[8] += deltaV[8];
		param[9] += deltaV[9];
		param[10] += deltaV[10];
		param[11] += deltaV[11];
		param[12] += deltaV[12];
		param[13] += deltaV[13];
		param[14] += deltaV[14];
		param[15] += deltaV[15];
		param[16] += deltaV[16];
		param[17] += deltaV[17];
		param[18] += deltaV[18];
		param[19] += deltaV[19];

		if (GetVecNorm(deltaV, 20) < eps)
		{
			break;
		}
	}

	delete[] deltaV;
}

template<class T>
void  __declspec(dllexport) Quad3dFitting(T param[18]
	, T(*point0)[3], T(*pointAff)[3], int point_num, double eps, int jt)
{
	T *deltaV = new T[18];
	T preNorm = INFINITYNUM;
	T rate = 1.0;
	for (int i = 0; i != jt; ++i)
	{
		GetDeltaVec4Quad3d(deltaV, param, point0, pointAff, point_num);

		param[0] += deltaV[0] * rate;
		param[1] += deltaV[1] * rate;
		param[2] += deltaV[2] * rate;
		param[3] += deltaV[3] * rate;
		param[4] += deltaV[4] * rate;
		param[5] += deltaV[5] * rate;
		param[6] += deltaV[6] * rate;
		param[7] += deltaV[7] * rate;
		param[8] += deltaV[8] * rate;
		param[9] += deltaV[9] * rate;
		param[10] += deltaV[10] * rate;
		param[11] += deltaV[11] * rate;
		param[12] += deltaV[12] * rate;
		param[13] += deltaV[13] * rate;
		param[14] += deltaV[14] * rate;
		param[15] += deltaV[15] * rate;
		param[16] += deltaV[16] * rate;
		param[17] += deltaV[17] * rate;
		T norm = GetVecNorm(deltaV, 18);

		if (norm < eps || norm>preNorm)
		{
			break;
		}

		preNorm = norm;
		rate /= 5;
	}

	delete[] deltaV;
}

template<class T>
void  __declspec(dllexport) QuadraticFitting(T param[12]
	, T(*point0)[3], T(*pointAff)[3], int point_num, double eps, int jt)
{
	T *deltaV = new T[12];
	for (int i = 0; i != jt; ++i)
	{
		getDeltaVec4Quadratic(deltaV, param, point0, pointAff, point_num);

		param[0] += deltaV[0];
		param[1] += deltaV[1];
		param[2] += deltaV[2];
		param[3] += deltaV[3];
		param[4] += deltaV[4];
		param[5] += deltaV[5];
		param[6] += deltaV[6];
		param[7] += deltaV[7];
		param[8] += deltaV[8];
		param[9] += deltaV[9];
		param[10] += deltaV[10];
		param[11] += deltaV[11];

		if (GetVecNorm(deltaV, 12) < eps)
		{
			break;
		}
	}

	delete[] deltaV;
}


template<class T>
void  __declspec(dllexport) PerspectMatFitting(T matP[9]
	, T(*point0)[3], T(*point)[3], int point_num, double eps, int jt)
{
	T *deltaV = new T[2];

	for (int i = 0; i != jt; ++i)
	{
		getDeltaVec4Perspect(deltaV, matP, point0, point, point_num);

		matP[6] += deltaV[0];
		matP[7] += deltaV[1];

		if (GetVecNorm(deltaV, 2) < eps)
		{
			break;
		}
	}

	delete[] deltaV;
}

template<class T>
void __declspec(dllexport) GetPlaneParam(T(&param)[3], T(*ps)[3], int num)
{
	planeFitting(param, ps, num, 1, 4);

	printf("plane's params is %f\t%f\t%f\n", param[0], param[1], param[2]);
}

template<class T>
void __declspec(dllexport) GetAffineTransformParam(T(&matrix)[9], T(*ps0)[3], T(*ps)[3], int num)
{
	T param[6] = { 0 };
	affineFitting(param, ps0, ps, num, 1, 4);

	//printf("%f\t%f\t%f\t%f\t%f\t%f\n", param[0], param[1], param[2], param[3], param[4], param[5]);
	memset(matrix, 0, sizeof(T) * 9);
	matrix[0] = param[0];
	matrix[1] = param[1];
	matrix[2] = param[4];
	matrix[3] = param[2];
	matrix[4] = param[3];
	matrix[5] = param[5];
	matrix[8] = param[1];
}


template<class T>
void __declspec(dllexport) GetAffine3dParamQR(T(&mat)[9], T(&t_a)[3]
	, T(&ortMat)[9], T(&tVec)[3], T(*ps0)[3], T(*ps)[3], int num)
{
	T param[9] = { 0 };
	Affine3dFittingQR(param, ortMat, tVec, ps0, ps, num, 1, 3);

	mat[0] = ortMat[0] * param[3];
	mat[1] = ortMat[0] * param[4] + ortMat[1] * param[6];
	mat[2] = ortMat[0] * param[5] + ortMat[1] * param[7] + ortMat[2] * param[8];
	mat[3] = ortMat[3] * param[3];
	mat[4] = ortMat[3] * param[4] + ortMat[4] * param[6];
	mat[5] = ortMat[3] * param[5] + ortMat[4] * param[7] + ortMat[5] * param[8];
	mat[6] = ortMat[6] * param[3];
	mat[7] = ortMat[6] * param[4] + ortMat[7] * param[6];
	mat[8] = ortMat[6] * param[5] + ortMat[7] * param[7] + ortMat[8] * param[8];

	t_a[0] = tVec[0] + param[0];
	t_a[1] = tVec[1] + param[1];
	t_a[2] = tVec[2] + param[2];
}



template<class T>
void __declspec(dllexport) GetAffineTransformParamQR(T(&matrix)[9], T(&ortMat)[9]
	, T(*ps0)[3], T(*ps)[3], int num)
{
	T param[5] = { 0 };
	AffineFittingQR(param, ortMat, ps0, ps, num, 1, 3);

	//printf("%f\t%f\t%f\t%f\t%f\t%f\n", param[0], param[1], param[2], param[3], param[4], param[5]);
	memset(matrix, 0, sizeof(T) * 9);
	matrix[0] = ortMat[0] * param[2];
	matrix[1] = ortMat[0] * param[3] + ortMat[1] * param[4];
	matrix[2] = param[0] + ortMat[2];
	matrix[3] = ortMat[3] * param[2];
	matrix[4] = ortMat[3] * param[3] + ortMat[4] * param[4];
	matrix[5] = param[1] + ortMat[5];
	matrix[8] = 1;
}

template<class T>
void __declspec(dllexport) GetAnalyTransParam(T(*mat_set)[MAXMATSIZE], int deg
	, T(*ps0)[3], T(*ps)[3], int num)
{
	//ADMMFitting(mat_set, deg, ps0, ps, num, 0.01, 50);
	
	for (int i = 1; i <= deg; ++i)
	{
		AnalyticFitting(mat_set, i, ps0, ps, num, 0.01, 3);
	}
}

template<class T>
void __declspec(dllexport) GetCubicTransformParam(T(&matC)[9], T(&matQ)[9], T(&linearT)[9]
	, T(*ps0)[3], T(*ps)[3], int num)
{
	T param[20] = { 0 };
	cubicFitting(param, ps0, ps, num, 0.01, 50);

	memset(matC, 0, sizeof(T) * 9);
	matC[0] = param[12];
	matC[1] = param[13];
	matC[2] = param[14];
	matC[3] = param[15];
	matC[4] = param[16];
	matC[5] = param[17];
	matC[6] = param[18];
	matC[7] = param[19];

	memset(matQ, 0, sizeof(T) * 9);
	matQ[0] = param[6];
	matQ[1] = param[7];
	matQ[2] = param[8];
	matQ[3] = param[9];
	matQ[4] = param[10];
	matQ[5] = param[11];

	memset(linearT, 0, sizeof(T) * 9);
	linearT[0] = param[2];
	linearT[1] = param[3];
	linearT[2] = param[0];
	linearT[3] = param[4];
	linearT[4] = param[5];
	linearT[5] = param[1];
	linearT[8] = 1;
}

template<class T>
void __declspec(dllexport) GetQuadTransf3dParam(T(&coef2)[18]
	, T(*ps0)[3], T(*ps)[3], int num)
{
	T param[18] = { 0 };
	Quad3dFitting(param, ps0, ps, num, 0.00001, 64);

	coef2[0] = param[0];
	coef2[1] = param[3];
	coef2[2] = param[6];
	coef2[3] = param[9];
	coef2[4] = param[12];
	coef2[5] = param[15];

	coef2[6] = param[1];
	coef2[7] = param[4];
	coef2[8] = param[7];
	coef2[9] = param[10];
	coef2[10] = param[13];
	coef2[11] = param[16];

	coef2[12] = param[2];
	coef2[13] = param[5];
	coef2[14] = param[8];
	coef2[15] = param[11];
	coef2[16] = param[14];
	coef2[17] = param[17];
}

template<class T>
void __declspec(dllexport) GetQuadraticTransformParam(T(&matQ)[9], T(&linearT)[9]
	, T(*ps0)[3], T(*ps)[3], int num)
{
	T param[12] = { 0 };
	quadraticFitting(param, ps0, ps, num, 0.01, 50);

	memset(matQ, 0, sizeof(T) * 9);
	matQ[0] = param[6];
	matQ[1] = param[7];
	matQ[2] = param[8];
	matQ[3] = param[9];
	matQ[4] = param[10];
	matQ[5] = param[11];

	memset(linearT, 0, sizeof(T) * 9);
	linearT[0] = param[2];
	linearT[1] = param[3];
	linearT[2] = param[0];
	linearT[3] = param[4];
	linearT[4] = param[5];
	linearT[5] = param[1];
	linearT[8] = 1;
}

template<class T>
void __declspec(dllexport) GetPerspectiveTransform(T(&matP)[9]
	, T(*ps0)[3], T(*ps)[3], int num)
{
	perspectMatFitting(matP, ps0, ps, num, 0.001, 5);

	//printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", matP[0], matP[1], matP[2], matP[3]
	//	, matP[4], matP[5], matP[6], matP[7], matP[8]);
}


inline double Analytic_ICP_3D(CPoint3 *point_set4regist, int num4regist
	, CPoint3 *aim_point_set, int aim_num, int icpIterCount, double eps)
{
	vector<int> fit_set;

	CPoint3 center, aimCenter;
	double dist;
	dist = 0;
	for (int i = 0; i != num4regist; ++i)
	{
		dist += sqrt(pow(aim_point_set[i][0] - point_set4regist[i][0], 2)
			+ pow(aim_point_set[i][1] - point_set4regist[i][1], 2)
			+ pow(aim_point_set[i][2] - point_set4regist[i][2], 2));

	}
	dist /= num4regist;

	printf("initial error: %.6f\n\n", dist);
	for (int i = 0; i != aim_num; ++i)
	{
		aimCenter[0] += aim_point_set[i][0] / aim_num;
		aimCenter[1] += aim_point_set[i][1] / aim_num;
		aimCenter[2] += aim_point_set[i][2] / aim_num;
	}

	for (int i = 0; i != num4regist; ++i)
	{
		center[0] += point_set4regist[i][0] / num4regist;
		center[1] += point_set4regist[i][1] / num4regist;
		center[2] += point_set4regist[i][2] / num4regist;
	}
	CPoint3 trans = aimCenter - center;
	for (int i = 0; i != num4regist; ++i)
	{
		point_set4regist[i][0] += trans[0];
		point_set4regist[i][1] += trans[1];
		point_set4regist[i][2] += trans[2];
	}

	dist = 0;
	for (int i = 0; i != num4regist; ++i)
	{
		dist += sqrt(pow(aim_point_set[i][0] - point_set4regist[i][0], 2)
			+ pow(aim_point_set[i][1] - point_set4regist[i][1], 2)
			+ pow(aim_point_set[i][2] - point_set4regist[i][2], 2));

	}
	dist /= num4regist;

	printf("initial after error: %.6f\n\n", dist);

	Matrix R = Matrix::eye(3);
	Matrix t(3, 1);
	double rigidErr;
	
	double preDist = INFINITYNUM;
	double result = 0;
	double icpICount = icpIterCount;
	IcpPointToPoint icp(aim_point_set, aim_num, icpIterCount);

	double18 coef2;
	double9 coef1;
	double3 coef0;
	memset(coef2, 0, sizeof(double18));
	memset(coef1, 0, sizeof(double9));
	memset(coef0, 0, sizeof(double3));

	double3 *aim_d_set = new double3[num4regist];
	double3 *d_set = new double3[num4regist];
	double3 d_t;

	FILE *fp = fopen("F:\\analytic3d_error_time.csv", "w+");
	if (fp == NULL) {
		fprintf(stderr, "fopen() failed.\n");
		exit(EXIT_FAILURE);
	}

	clock_t start, finish, icpStart, icpFinish, finishStep;

	double duration, durationStep, icpDuration;

	start = clock();

	for (int i = 0; i != MAXNEWICPITERCOUNT; ++i)
	{
		icpStart = clock();
		icp.setMaxIterations(icpICount);
		icp.fit(point_set4regist, num4regist, R, t, fit_set, -1);

		for (int j = 0; j != fit_set.size(); ++j)
		{
			matrix3double3(aim_d_set[j], aim_point_set[fit_set[j]]);
		}

		if (0 == i)
		{
			icpFinish = clock();
			icpDuration = (double)(icpFinish - icpStart) / CLOCKS_PER_SEC;

			printf("icp running time: %.5f s\n", icpDuration);
		}

		//for (int j = 0; j != num4regist; ++j)
		//{
		//	//momentum
		//	point_set4regist[j] = R*point_set4regist[j] + t;
		//}

		matrix3double3(d_set, point_set4regist, num4regist);

		double9 ort;
		double3 tV;

		ort[0] = R.val[0][0];
		ort[1] = R.val[0][1];
		ort[2] = R.val[0][2];
		ort[3] = R.val[1][0];
		ort[4] = R.val[1][1];
		ort[5] = R.val[1][2];
		ort[6] = R.val[2][0];
		ort[7] = R.val[2][1];
		ort[8] = R.val[2][2];

		tV[0] = t.val[0][0];
		tV[1] = t.val[1][0];
		tV[2] = t.val[2][0];
		GetAffine3dParamQR(coef1, coef0, ort, tV, d_set, aim_d_set, num4regist);

		coef0[0] += tV[0];
		coef0[1] += tV[1];
		coef0[2] += tV[2];

		coefmat *mat_set = new coefmat[MAXDEGREE];
		memset(mat_set, 0, sizeof(coefmat)*MAXDEGREE);
		mat_set[0][0] = coef0[0];
		mat_set[0][1] = coef0[1];
		mat_set[0][2] = coef0[2];

		mat_set[1][0] = coef1[0];
		mat_set[1][1] = coef1[1];
		mat_set[1][2] = coef1[2];
		mat_set[1][3] = coef1[3];
		mat_set[1][4] = coef1[4];
		mat_set[1][5] = coef1[5];
		mat_set[1][6] = coef1[6];
		mat_set[1][7] = coef1[7];
		mat_set[1][8] = coef1[8];

		dist = 0;
		for (int i = 0; i != num4regist; ++i)
		{
			AnalyticTrans3d_1(mat_set, d_set[i], d_t);
			double2matrix(point_set4regist[i], d_t);

			aim_d_set[i][0] -= d_t[0];
			aim_d_set[i][1] -= d_t[1];
			aim_d_set[i][2] -= d_t[2];

			dist += sqrt(pow(aim_d_set[i][0], 2) 
				+ pow(aim_d_set[i][1], 2) 
				+ pow(aim_d_set[i][2], 2));

		}
		dist /= num4regist;

		if (i > 1)
		{
			GetQuadTransf3dParam(coef2, d_set, aim_d_set, num4regist);

			mat_set[2][0] = coef2[0];
			mat_set[2][1] = coef2[1];
			mat_set[2][2] = coef2[2];
			mat_set[2][3] = coef2[3];
			mat_set[2][4] = coef2[4];
			mat_set[2][5] = coef2[5];
			mat_set[2][6] = coef2[6];
			mat_set[2][7] = coef2[7];
			mat_set[2][8] = coef2[8];
			mat_set[2][9] = coef2[9];
			mat_set[2][10] = coef2[10];
			mat_set[2][11] = coef2[11];
			mat_set[2][12] = coef2[12];
			mat_set[2][13] = coef2[13];
			mat_set[2][14] = coef2[14];
			mat_set[2][15] = coef2[15];
			mat_set[2][16] = coef2[16];
			mat_set[2][17] = coef2[17];
			dist = 0;

			for (int i = 0; i != num4regist; ++i)
			{
				AnalyticTrans3d_2(mat_set[2], d_set[i], d_t);
				point_set4regist[i][0] += d_t[0];
				point_set4regist[i][1] += d_t[1];
				point_set4regist[i][2] += d_t[2];

				dist += get_euclid_dist<double, 3>(d_t, aim_d_set[i]);
			}

			dist /= num4regist;
		}

	

		if (icpICount > 2)
		{
			icpICount /= 1.5;
		}

		if (dist > preDist || dist < eps || dist < 0.02)
		{
			result = dist;

			finishStep = clock();
			durationStep = (double)(finishStep - start) / CLOCKS_PER_SEC;
			fprintf(fp, "%f,%f\n", (float)dist, durationStep);

			printf("-------------final error is %.5f-------------\n", (float)dist);
			break;
		}

		finishStep = clock();
		durationStep = (double)(finishStep - start) / CLOCKS_PER_SEC;
		fprintf(fp, "%f,%f\n", (float)dist, durationStep);
		
		preDist = dist;


	}
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf("ftransform algorithms' response time is--------------%.5f seconds-------------", duration);

	fclose(fp);
	printf("%f\t%f\t%f\n", coef0[0], coef0[1], coef0[2]);
	printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"
		, coef1[0], coef1[1], coef1[2]
		, coef1[3], coef1[4], coef1[5]
		, coef1[6], coef1[7], coef1[8]);

	printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"
		, coef2[0], coef2[1], coef2[2]
		, coef2[3], coef2[4], coef2[5]
		, coef2[6], coef2[7], coef2[8]
		, coef2[9], coef2[10], coef2[11]
		, coef2[12], coef2[13], coef2[14]
		, coef2[15], coef2[16], coef2[17]);
	delete[]aim_d_set;
	delete[] d_set;
	return result;
}


/*!
* Non-rigid registration for Small nonlinear deformation.
* \param point_set4regist
* \param num4regist
* \param aim_point_set
* \param aim_num
* \param icpIterCount
* \param eps accuracy
* \return score for several-for-one
*/
inline double Analytic_ICP_2D(CPoint2 *point_set4regist, int num4regist
	, CPoint2 *aim_point_set, int aim_num, int icpIterCount, double eps)
{
	vector<int> fit_set;

	CPoint2 center, aimCenter;
	for (int i = 0; i != aim_num; ++i)
	{
		aimCenter[0] += aim_point_set[i][0] / aim_num;
		aimCenter[1] += aim_point_set[i][1] / aim_num;
	}

	for (int i = 0; i != num4regist; ++i)
	{
		center[0] += point_set4regist[i][0] / num4regist;
		center[1] += point_set4regist[i][1] / num4regist;
	}
	CPoint2 trans = aimCenter - center;
	for (int i = 0; i != num4regist; ++i)
	{
		point_set4regist[i][0] += trans[0];
		point_set4regist[i][1] += trans[1];
	}

	Matrix R = Matrix::eye(2);
	Matrix t(2, 1);
	double rigidErr;
	double dist;
	double preDist = INFINITYNUM;
	double result = 0;
	double icpICount = icpIterCount;
	IcpPointToPoint icp(aim_point_set, aim_num, icpIterCount);
	double9 matrix;
	double9 matQ;
	double9 matC;
	double3 *aim_d_set = new double3[num4regist];
	double3 *d_set = new double3[num4regist];
	double3 d_t;

	FILE *fp = fopen("F:\\icp_error_time.csv", "w+");
	if (fp == NULL) {
		fprintf(stderr, "fopen() failed.\n");
		exit(EXIT_FAILURE);
	}

	clock_t start, finish, icpStart, icpFinish,finishStep;

	double duration, durationStep, icpDuration;

	start = clock();

	for (int i = 0; i != MAXNEWICPITERCOUNT; ++i)
	{
		icpStart = clock();
		icp.setMaxIterations(icpICount);
		icp.fit(point_set4regist, num4regist, R, t, fit_set, -1);

		if (0 == i)
		{
			rigidErr = 0;
			for (int j = 0; j != num4regist; ++j)
			{
				matrix2double3(d_t, point_set4regist[j]);
				rigidErr += get_euclid_dist<double, 2>(d_t, aim_d_set[j]);
			}

			rigidErr /= num4regist;

			printf("-------------rigid error is %f-------------\n", (float)rigidErr);
		}
		for (int j = 0; j != num4regist; ++j)
		{
			point_set4regist[j] = R*point_set4regist[j] + t;
		}

		if (0 == i)
		{
			icpFinish = clock();
			icpDuration = (double)(icpFinish - icpStart) / CLOCKS_PER_SEC;

			printf("icp running time: %.5f s\n", icpDuration);
		}

		for (int j = 0; j != fit_set.size(); ++j)
		{
			matrix2double3(aim_d_set[j], aim_point_set[fit_set[j]]);
		}


		for (int j = 0; j != num4regist; ++j)
		{
			//momentum
			point_set4regist[j] = R*point_set4regist[j] + t;
		}

		matrix2double3(d_set, point_set4regist, num4regist);

		double9 ort;
		memset(ort, 0, sizeof(double9));
		ort[0] = R.val[0][0];
		ort[1] = R.val[0][1];
		ort[2] = t.val[0][0];
		ort[3] = R.val[1][0];
		ort[4] = R.val[1][1];
		ort[5] = t.val[1][0];
		ort[8] = 1;
		GetAffineTransformParamQR(matrix, ort, d_set, aim_d_set, num4regist);
		//GetPerspectiveTransform(matrix, d_set, aim_d_set, num4regist);
		//GetQuadraticTransformParam(matQ, matrix, d_set, aim_d_set, num4regist);
		//GetCubicTransformParam(matC, matQ, matrix, d_set, aim_d_set, num4regist);
		//memset(matQ, 0, sizeof(double9));
		dist = 0;

		coefmat *mat_set = new coefmat[MAXDEGREE];
		memset(mat_set, 0, sizeof(coefmat)*MAXDEGREE);
		mat_set[0][0] = matrix[2];
		mat_set[0][1] = matrix[5];

		mat_set[1][0] = matrix[0];
		mat_set[1][1] = matrix[1];
		mat_set[1][2] = matrix[3];
		mat_set[1][3] = matrix[4];

		GetAnalyTransParam(mat_set, 3, d_set, aim_d_set, num4regist);

		for (int i = 0; i != num4regist; ++i)
		{
			AnalyticTrans(mat_set, 3, d_set[i], d_t);
			double2matrix(point_set4regist[i], d_t);

			dist += get_euclid_dist<double, 2>(d_t, aim_d_set[i]);
		}

		dist /= num4regist;

		if (icpICount > 2)
		{
			icpICount /= 1.5;
		}

		if (dist > preDist && dist < eps)
		{
			result = dist;

			printf("-------------final error is %.5f-------------\n", (float)dist);
			break;
		}

		preDist = dist;

		finishStep = clock();
		durationStep = (double)(finishStep - start) / CLOCKS_PER_SEC;
		fprintf(fp, "%f,%f\n", (float)dist, durationStep);
	}
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf("ftransform algorithms' response time is--------------%.5f seconds-------------", duration);

	fclose(fp);
	printf("%f\t%f\t%f\t%f\t%f\t%f\n"
		, matrix[0], matrix[1], matrix[2], matrix[3], matrix[4], matrix[5]);
	delete[]aim_d_set;
	delete[] d_set;
	return result;
}


#endif
