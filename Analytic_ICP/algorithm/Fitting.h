#ifndef FITTING
#define FITTING
/*!
*      \file Fitting.h
*      \brief Fitting
*	   \author Wei Feng
*      \date 10/08/2020
*
*/
#include "Algo.h"

#define MAXNEWICPITERCOUNT 14

template<class T>
void  __declspec(dllexport) GetMatB4Affine(MatrixXd &matB, T(*point0)[3], T(*point)[3], int point_num)
{
	//memset(matB, 0, sizeof(T) * 2 * 6 * point_num);

	for (int i = 0; i != point_num; ++i)
	{
		matB(i * 2 + 0, 0) = point0[i][0];
		matB(i * 2 + 0, 1) = point0[i][1];
		matB(i * 2 + 0, 2) = 0;
		matB(i * 2 + 0, 3) = 0;
		matB(i * 2 + 0, 4) = 1;
		matB(i * 2 + 0, 5) = 0;

		matB(i * 2 + 1, 0) = 0;
		matB(i * 2 + 1, 1) = 0;
		matB(i * 2 + 1, 2) = point0[i][0];
		matB(i * 2 + 1, 3) = point0[i][1];
		matB(i * 2 + 1, 4) = 0;
		matB(i * 2 + 1, 5) = 1;
	}
}

template<class T>
void  __declspec(dllexport) GetMatB3d4AffineQR(MatrixXd &matB, T(&ortMat)[9]
	, T(*point0)[3], int point_num)
{
	//memset(matB, 0, sizeof(T) * 3 * 9 * point_num);

	for (int i = 0; i != point_num; ++i)
	{
		matB(i * 3, 0) = 1;
		matB(i * 3, 1) = 0;
		matB(i * 3, 2) = 0;
		matB(i * 3, 3) = ortMat[0] * point0[i][0];
		matB(i * 3, 4) = ortMat[0] * point0[i][1];
		matB(i * 3, 5) = ortMat[0] * point0[i][2];
		matB(i * 3, 6) = ortMat[1] * point0[i][1];
		matB(i * 3, 7) = ortMat[1] * point0[i][2];
		matB(i * 3, 8) = ortMat[2] * point0[i][2];

		matB(i * 3 + 1, 0) = 0;
		matB(i * 3 + 1, 1) = 1;
		matB(i * 3 + 1, 2) = 0;
		matB(i * 3 + 1, 3) = ortMat[3] * point0[i][0];
		matB(i * 3 + 1, 4) = ortMat[3] * point0[i][1];
		matB(i * 3 + 1, 5) = ortMat[3] * point0[i][2];
		matB(i * 3 + 1, 6) = ortMat[4] * point0[i][1];
		matB(i * 3 + 1, 7) = ortMat[4] * point0[i][2];
		matB(i * 3 + 1, 8) = ortMat[5] * point0[i][2];

		matB(i * 3 + 2, 0) = 0;
		matB(i * 3 + 2, 1) = 0;
		matB(i * 3 + 2, 2) = 1;
		matB(i * 3 + 2, 3) = ortMat[6] * point0[i][0];
		matB(i * 3 + 2, 4) = ortMat[6] * point0[i][1];
		matB(i * 3 + 2, 5) = ortMat[6] * point0[i][2];
		matB(i * 3 + 2, 6) = ortMat[7] * point0[i][1];
		matB(i * 3 + 2, 7) = ortMat[7] * point0[i][2];
		matB(i * 3 + 2, 8) = ortMat[8] * point0[i][2];
	}
}



template<class T>
void  __declspec(dllexport) GetMatB4AffineQR(MatrixXd &matB, T(&ortMat)[9]
	, T(*point0)[3], int point_num)
{
	//memset(matB, 0, sizeof(T) * 2 * 5 * point_num);

	for (int i = 0; i != point_num; ++i)
	{
		matB(i * 2, 0) = 1;
		matB(i * 2, 1) = 0;
		matB(i * 2, 2) = ortMat[0] * point0[i][0];
		matB(i * 2, 3) = ortMat[0] * point0[i][1];
		matB(i * 2, 4) = ortMat[1] * point0[i][1];

		matB(i * 2 + 1, 0) = 0;
		matB(i * 2 + 1, 1) = 1;
		matB(i * 2 + 1, 2) = ortMat[3] * point0[i][0];
		matB(i * 2 + 1, 3) = ortMat[3] * point0[i][1];
		matB(i * 2 + 1, 4) = ortMat[4] * point0[i][1];
	}
}


template<class T>
void  __declspec(dllexport) GetMatB4Nominal(MatrixXd &matB, int deg, T(*mat_set)[MAXMATSIZE],
	T nomiVal[3], T(*point0)[3], int point_num)
{
	int i, j, k;
	double one_fac;
	int comb_num;
	for (i = 0; i != point_num; ++i)
	{
		for (j = 1; j <= deg; ++j)
		{
			one_fac = 1.0 / Factorial(j);
			for (k = 0; k <= j; ++k)
			{
				comb_num = Combination(j, k);
				if (j - k >= 1)
				{
					matB(i * 2 + 0, 0) -= one_fac*comb_num*(j - k)
						*pow(point0[i][0] - nomiVal[0], j - k - 1)
						*pow(point0[i][1] - nomiVal[1], k)*mat_set[j][k];

					matB(i * 2 + 1, 0) -= one_fac*comb_num*(j - k)
						*pow(point0[i][0] - nomiVal[0], j - k - 1)
						*pow(point0[i][1] - nomiVal[1], k)*mat_set[j][j + k + 1];
				}

				if (k >= 1)
				{
					matB(i * 2 + 0, 1) -= one_fac*comb_num*k
						*pow(point0[i][0] - nomiVal[0], j - k)
						*pow(point0[i][1] - nomiVal[1], k - 1)*mat_set[j][k];

					matB(i * 2 + 1, 1) -= one_fac*comb_num*k
						*pow(point0[i][0] - nomiVal[0], j - k)
						*pow(point0[i][1] - nomiVal[1], k - 1)*mat_set[j][j + k + 1];
				}

			}
		}
	}
}


template<class T>
void  __declspec(dllexport) GetMatB4Analyt3d(MatrixXd &matB, int deg
	, T(*point0)[3], int point_num)
{
	int param_num = 0;
	int coef_size;
	double coef[MAXNOMIALSIZE] = { 0 };
	int i, j, k;
	for (i = 0; i != point_num; ++i)
	{
		param_num = 0;
		for (j = 0; j <= deg; ++j)
		{
			GetCoefOfTaylor(coef, coef_size, point0[i], j);
		
			for (int k = 0; k != coef_size; ++k)
			{
				matB(i * 3 + 0, param_num + k) = coef[k];
				matB(i * 3 + 0, param_num + coef_size + k) = 0;
				matB(i * 3 + 0, param_num + coef_size * 2 + k) = 0;


				matB(i * 3 + 1, param_num + k) = 0;
				matB(i * 3 + 1, param_num + coef_size + k) = coef[k];
				matB(i * 3 + 1, param_num + coef_size * 2 + k) = 0;

				matB(i * 3 + 2, param_num + k) = 0;
				matB(i * 3 + 2, param_num + coef_size + k) = 0;
				matB(i * 3 + 2, param_num + coef_size * 2 + k) = coef[k];
			}
			param_num += coef_size * 3;

		}
	}
}


template<class T>
void  __declspec(dllexport) GetMatB4Analytic(MatrixXd &matB, int deg, int c,
	T nomiVal[3], T(*point0)[3], int point_num)
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
				matB(i * 2 + 0, e_i + k) = one_fac*comb_num
					*pow(point0[i][0] - nomiVal[0], j - k)
					*pow(point0[i][1] - nomiVal[1], k);
				matB(i * 2 + 1, e_i + k) = 0;

				matB(i * 2 + 0, e_i + j + 1 + k) = 0;
				matB(i * 2 + 1, e_i + j + 1 + k) = one_fac*comb_num
					*pow(point0[i][0] - nomiVal[0], j - k)
					*pow(point0[i][1] - nomiVal[1], k);
			}
			e_i += 2 * (j + 1);
		}
	}
}


template<class T>
void  __declspec(dllexport) GetMatB4Quad3d(MatrixXd &matB, T(*point0)[3], int point_num)
{
	for (int i = 0; i != point_num; ++i)
	{
		matB(i * 3 + 0, 0) = 1;
		matB(i * 3 + 0, 1) = 0;
		matB(i * 3 + 0, 2) = 0;
		matB(i * 3 + 0, 3) = point0[i][0];
		matB(i * 3 + 0, 4) = point0[i][1];
		matB(i * 3 + 0, 5) = point0[i][2];
		matB(i * 3 + 0, 6) = 0;
		matB(i * 3 + 0, 7) = 0;
		matB(i * 3 + 0, 8) = 0;
		matB(i * 3 + 0, 9) = 0;
		matB(i * 3 + 0, 10) = 0;
		matB(i * 3 + 0, 11) = 0;
		matB(i * 3 + 0, 12) = 0.5*pow(point0[i][0], 2);
		matB(i * 3 + 0, 13) = point0[i][0] * point0[i][1];
		matB(i * 3 + 0, 14) = point0[i][0] * point0[i][2];
		matB(i * 3 + 0, 15) = 0.5*pow(point0[i][1], 2);
		matB(i * 3 + 0, 16) = point0[i][1] * point0[i][2];
		matB(i * 3 + 0, 17) = 0.5*pow(point0[i][2], 2);
		matB(i * 3 + 0, 18) = 0;
		matB(i * 3 + 0, 19) = 0;
		matB(i * 3 + 0, 20) = 0;
		matB(i * 3 + 0, 21) = 0;
		matB(i * 3 + 0, 22) = 0;
		matB(i * 3 + 0, 23) = 0;
		matB(i * 3 + 0, 24) = 0;
		matB(i * 3 + 0, 25) = 0;
		matB(i * 3 + 0, 26) = 0;
		matB(i * 3 + 0, 27) = 0;
		matB(i * 3 + 0, 28) = 0;
		matB(i * 3 + 0, 29) = 0;



		matB(i * 3 + 1, 0) = 0;
		matB(i * 3 + 1, 1) = 1;
		matB(i * 3 + 1, 2) = 0;
		matB(i * 3 + 1, 3) = 0;
		matB(i * 3 + 1, 4) = 0;
		matB(i * 3 + 1, 5) = 0;
		matB(i * 3 + 1, 6) = point0[i][0];
		matB(i * 3 + 1, 7) = point0[i][1];
		matB(i * 3 + 1, 8) = point0[i][2];
		matB(i * 3 + 1, 9) = 0;
		matB(i * 3 + 1, 10) = 0;
		matB(i * 3 + 1, 11) = 0;
		matB(i * 3 + 1, 12) = 0;
		matB(i * 3 + 1, 13) = 0;
		matB(i * 3 + 1, 14) = 0;
		matB(i * 3 + 1, 15) = 0;
		matB(i * 3 + 1, 16) = 0;
		matB(i * 3 + 1, 17) = 0;
		matB(i * 3 + 1, 18) = 0.5*pow(point0[i][0], 2);
		matB(i * 3 + 1, 19) = point0[i][0] * point0[i][1];
		matB(i * 3 + 1, 20) = point0[i][0] * point0[i][2];
		matB(i * 3 + 1, 21) = 0.5*pow(point0[i][1], 2);
		matB(i * 3 + 1, 22) = point0[i][1] * point0[i][2];
		matB(i * 3 + 1, 23) = 0.5*pow(point0[i][2], 2);
		matB(i * 3 + 1, 24) = 0;
		matB(i * 3 + 1, 25) = 0;
		matB(i * 3 + 1, 26) = 0;
		matB(i * 3 + 1, 27) = 0;
		matB(i * 3 + 1, 28) = 0;
		matB(i * 3 + 1, 29) = 0;



		matB(i * 3 + 2, 0) = 0;
		matB(i * 3 + 2, 1) = 0;
		matB(i * 3 + 2, 2) = 1;
		matB(i * 3 + 2, 3) = 0;
		matB(i * 3 + 2, 4) = 0;
		matB(i * 3 + 2, 5) = 0;
		matB(i * 3 + 2, 6) = 0;
		matB(i * 3 + 2, 7) = 0;
		matB(i * 3 + 2, 8) = 0;
		matB(i * 3 + 2, 9) = point0[i][0];
		matB(i * 3 + 2, 10) = point0[i][1];
		matB(i * 3 + 2, 11) = point0[i][2];
		matB(i * 3 + 2, 12) = 0;
		matB(i * 3 + 2, 13) = 0;
		matB(i * 3 + 2, 14) = 0;
		matB(i * 3 + 2, 15) = 0;
		matB(i * 3 + 2, 16) = 0;
		matB(i * 3 + 2, 17) = 0;
		matB(i * 3 + 2, 18) = 0;
		matB(i * 3 + 2, 19) = 0;
		matB(i * 3 + 2, 20) = 0;
		matB(i * 3 + 2, 21) = 0;
		matB(i * 3 + 2, 22) = 0;
		matB(i * 3 + 2, 23) = 0;
		matB(i * 3 + 2, 24) = 0.5*pow(point0[i][0], 2);
		matB(i * 3 + 2, 25) = point0[i][0] * point0[i][1];
		matB(i * 3 + 2, 26) = point0[i][0] * point0[i][2];
		matB(i * 3 + 2, 27) = 0.5*pow(point0[i][1], 2);
		matB(i * 3 + 2, 28) = point0[i][1] * point0[i][2];
		matB(i * 3 + 2, 29) = 0.5*pow(point0[i][2], 2);
		
	}
}



template<class T>
void  __declspec(dllexport) GetL4Affine(VectorXd &l, T param[], T(*point0)[3], T(*point)[3], int point_num)
{
	//memset(l, 0, sizeof(T)*point_num * 2);

	for (int i = 0; i != point_num; ++i)
	{
		l(2 * i + 0) = point[i][0] - param[0] * point0[i][0] - param[1] * point0[i][1] - param[4];
		l(2 * i + 1) = point[i][1] - param[2] * point0[i][0] - param[3] * point0[i][1] - param[5];
	}
}

template<class T>
void  __declspec(dllexport) GetL3d4AffineQR(VectorXd &l, T(&ortMat)[9], T(&tVec)[3]
	, T param[], T(*point0)[3], T(*point)[3], int point_num)
{
	//memset(l, 0, sizeof(T)*point_num * 3);

	for (int i = 0; i != point_num; ++i)
	{
		l(3 * i + 0) = point[i][0] - tVec[0] - param[0]
			- ortMat[0] * param[3] * point0[i][0]
			- ortMat[0] * param[4] * point0[i][1]
			- ortMat[1] * param[6] * point0[i][1]
			- ortMat[0] * param[5] * point0[i][2]
			- ortMat[1] * param[7] * point0[i][2] 
			- ortMat[2] * param[8] * point0[i][2];
		l(3 * i + 1) = point[i][1] - tVec[1] - param[1]
			- ortMat[3] * param[3] * point0[i][0]
			- ortMat[3] * param[4] * point0[i][1]
			- ortMat[4] * param[6] * point0[i][1]
			- ortMat[3] * param[5] * point0[i][2]
			- ortMat[4] * param[7] * point0[i][2]
			- ortMat[5] * param[8] * point0[i][2];
		l(3 * i + 2) = point[i][2] - tVec[2] - param[2]
			- ortMat[6] * param[3] * point0[i][0]
			- ortMat[6] * param[4] * point0[i][1]
			- ortMat[7] * param[6] * point0[i][1]
			- ortMat[6] * param[5] * point0[i][2]
			- ortMat[7] * param[7] * point0[i][2]
			- ortMat[8] * param[8] * point0[i][2];
	}
}


template<class T>
void  __declspec(dllexport) GetL4AffineQR(VectorXd &l, T(&ortMat)[9]
	, T param[], T(*point0)[3], T(*point)[3], int point_num)
{
	//memset(l, 0, sizeof(T)*point_num * 2);

	for (int i = 0; i != point_num; ++i)
	{
		l(2 * i + 0) = point[i][0] - ortMat[2] - param[0]
			- ortMat[0] * param[2] * point0[i][0]
			- ortMat[0] * param[3] * point0[i][1]
			- ortMat[1] * param[4] * point0[i][1];
		l(2 * i + 1) = point[i][1] - ortMat[5] - param[1]
			- ortMat[3] * param[2] * point0[i][0]
			- ortMat[3] * param[3] * point0[i][1]
			- ortMat[4] * param[4] * point0[i][1];
	}
}

template<class T>
void  __declspec(dllexport) GetL4Analyt3d(VectorXd &l, VectorXd &param
	, int deg, T(*point0)[3], T(*point)[3], int point_num)
{
	int i, j, k;
	double coef[MAXNOMIALSIZE] = { 0 };
	int coef_size = 0;
	int temp_sn;
	for (i = 0; i != point_num; ++i)
	{
		l(3 * i + 0) = point[i][0];
		l(3 * i + 1) = point[i][1];
		l(3 * i + 2) = point[i][2];

		temp_sn = 0;

		for (j = 0; j <= deg; ++j)
		{
			GetCoefOfTaylor(coef, coef_size, point0[i], j);

			for (int k = 0; k != coef_size; ++k)
			{
				l(3 * i + 0) -= coef[k] * param(temp_sn + k);
				l(3 * i + 1) -= coef[k] * param(temp_sn + coef_size + k);
				l(3 * i + 2) -= coef[k] * param(temp_sn + coef_size * 2 + k);
			}

			temp_sn += coef_size * 3;
		}
	}
}

template<class T>
void  __declspec(dllexport) GetL4Analytic(VectorXd &l, T(*mat_set)[MAXMATSIZE], int deg
	, T nomiVal[3], T(*point0)[3], T(*point)[3], int point_num)
{
	//memset(l, 0, sizeof(T)*point_num * 2);
	int i, j, k;
	double one_fac;
	int comb_num;
	for (i = 0; i != point_num; ++i)
	{
		l(2 * i + 0) = point[i][0];
		l(2 * i + 1) = point[i][1];

		for (j = 0; j <= deg; ++j)
		{
			one_fac = 1.0 / Factorial(j);
			for (k = 0; k <= j; ++k)
			{
				comb_num = Combination(j, k);
				l(2 * i + 0) -= one_fac*comb_num
					*pow(point0[i][0] - nomiVal[0], j - k)
					*pow(point0[i][1] - nomiVal[1], k)*mat_set[j][k];
				l(2 * i + 1) -= one_fac*comb_num
					*pow(point0[i][0] - nomiVal[0], j - k)
					*pow(point0[i][1] - nomiVal[1], k)*mat_set[j][j + k + 1];
			}
		}
	}
}


template<class T>
void  __declspec(dllexport) GetL4Quad3d(VectorXd &l
	, VectorXd &param, T(*point0)[3], T(*point)[3], int point_num)
{
	for (int i = 0; i != point_num; ++i)
	{
		l(3 * i + 0) = point[i][0] - param(0) 
			- param(3)*point0[i][0] 
			- param(4)*point0[i][1] 
			- param(5)*point0[i][2]
			-0.5*param(12) * pow(point0[i][0], 2)
			- param(13) * point0[i][0] * point0[i][1]
			- param(14) * point0[i][0] * point0[i][2]
			- 0.5*param(15) * pow(point0[i][1], 2)
			- param(16) * point0[i][1] * point0[i][2]
			- 0.5*param(17) * pow(point0[i][2], 2);

		l(3 * i + 1) = point[i][1] - param(1)
			- param(6)*point0[i][0]
			- param(7)*point0[i][1]
			- param(8)*point0[i][2]
			- 0.5*param(18) * pow(point0[i][0], 2)
			- param(19) * point0[i][0] * point0[i][1]
			- param(20) * point0[i][0] * point0[i][2]
			- 0.5*param(21) * pow(point0[i][1], 2)
			- param(22) * point0[i][1] * point0[i][2]
			- 0.5*param(23) * pow(point0[i][2], 2);

		l(3 * i + 2) = point[i][2] - param(2)
			- param(9)*point0[i][0]
			- param(10)*point0[i][1]
			- param(11)*point0[i][2]
			- 0.5*param(24) * pow(point0[i][0], 2)
			- param(25) * point0[i][0] * point0[i][1]
			- param(26) * point0[i][0] * point0[i][2]
			- 0.5*param(27) * pow(point0[i][1], 2)
			- param(28) * point0[i][1] * point0[i][2]
			- 0.5*param(29) * pow(point0[i][2], 2);
	}
}


template<class T>
void  __declspec(dllexport) GetDeltaVec4Affine(VectorXd &DeltaV, T param[]
	, T(*point0)[3], T(*point)[3], int point_num)
{
	Matrix2d matB = MatrixXd::Zero(2 * point_num, 6);
	VectorXd l = VectorXd::Zero(point_num * 2);

	GetMatB4Affine(matB, point0, point, point_num);
	GetL4Affine(l, param, point0, point, point_num);
	GetCorrection(DeltaV, matB, l, point_num * 2, 6);
}


template<class T>
void  __declspec(dllexport) GetDVec3d4AffineQR(VectorXd &DVec, T(&ortMat)[9], T(&tVec)[3]
	, T param[], T(*point0)[3], T(*point)[3], int point_num)
{
	MatrixXd matB = MatrixXd::Zero(3 * point_num, 9);
	VectorXd l = VectorXd::Zero(point_num * 3);

	GetMatB3d4AffineQR(matB, ortMat, point0, point_num);
	GetL3d4AffineQR(l, ortMat, tVec, param, point0, point, point_num);
	GetCorrection(DVec, matB, l, point_num * 3, 9);
}

template<class T>
void  __declspec(dllexport) GetDeltaVec4AffineQR(VectorXd &DeltaV, T(&ortMat)[9], T param[]
	, T(*point0)[3], T(*point)[3], int point_num)
{
	MatrixXd matB = MatrixXd::Zero(2 * point_num, 5);
	VectorXd l = VectorXd::Zero(point_num * 2);

	GetMatB4AffineQR(matB, ortMat, point0, point_num);
	GetL4AffineQR(l, ortMat, param, point0, point, point_num);
	GetCorrection(DeltaV, matB, l, point_num * 2, 5);
}


template<class T>
T  __declspec(dllexport) GetDeltaVec4Nominal(VectorXd &deltaV, T(*mat_set)[MAXMATSIZE]
	, T nomiVal[3]
	, int deg, T(*point0)[3], T(*point)[3], int point_num)
{
#define D 2
	MatrixXd matB = MatrixXd::Zero(2 * point_num, D);
	VectorXd l = VectorXd::Zero(point_num * 2);


	GetMatB4Nominal(matB, deg, mat_set, nomiVal, point0, point_num);
	GetL4Analytic(l, mat_set, deg, nomiVal, point0, point, point_num);
	GetCorrection(deltaV, matB, l, point_num * 2, D);

	T error = deltaV.norm();

	return error;
}


template<class T>
T  __declspec(dllexport) GetDeltaVec3d4Analytic(VectorXd &deltaV
	, int param_num, VectorXd &param
	, int deg, T(*point0)[3], T(*point)[3], int point_num)
{

	MatrixXd matB = MatrixXd::Zero(3 * point_num, param_num);
	VectorXd l = VectorXd::Zero(point_num * 3);

	GetMatB4Analyt3d(matB, deg, point0, point_num);
	GetL4Analyt3d(l, param, deg, point0, point, point_num);
	GetCorrection(deltaV, matB, l, point_num * 3, param_num);

	T error = deltaV.norm();

	return error;
}


template<class T>
T  __declspec(dllexport) GetDeltaVec4Analytic(T(*delta_set)[MAXMATSIZE], T(*mat_set)[MAXMATSIZE]
	, T nomiVal[3]
	, int deg, T(*point0)[3], T(*point)[3], int point_num)
{
	int sn = Sn(2, 2, deg + 1);
	MatrixXd matB = MatrixXd::Zero(2 * point_num, sn);
	VectorXd l = VectorXd::Zero(point_num * 2);
	VectorXd deltaV = VectorXd::Zero(sn);

	GetMatB4Analytic(matB, deg, sn, nomiVal, point0, point_num);
	GetL4Analytic(l, mat_set, deg, nomiVal, point0, point, point_num);
	GetCorrection(deltaV, matB, l, point_num * 2, sn);

	int sn_temp;
	for (int i = 0; i <= deg; ++i)
	{
		sn_temp = Sn(2, 2, i);

		for (int j = 0; j != 2 * (i + 1); ++j)
		{
			delta_set[i][j] = deltaV(sn_temp + j);
		}
	}

	T error = deltaV.norm();

	return error;
}

template<class T>
void  __declspec(dllexport) GetDeltaVec4Quad3d(VectorXd &DeltaV
	, MatrixXd &matB, VectorXd &l, VectorXd &param
	, T(*point0)[3], T(*pointAff)[3], int point_num)
{
	GetMatB4Quad3d(matB, point0, point_num);
	GetL4Quad3d(l, param, point0, pointAff, point_num);
	GetCorrection(DeltaV, matB, l, point_num * 3, 30);
	
}


template<class T>
void  __declspec(dllexport) AffineFitting(T param[6], T(*point0)[3], T(*point)[3]
	, int point_num, double eps, int jt)
{
	VectorXd deltaV = VectorXd::Zero(6);

	param[0] = 1;
	param[1] = 0;
	param[2] = 0;
	param[3] = 1;
	param[4] = 0;
	param[5] = 0;

	for (int i = 0; i != jt; ++i)
	{
		GetDeltaVec4Affine(deltaV, param, point0, point, point_num);

		param[0] += deltaV(0);
		param[1] += deltaV(1);
		param[2] += deltaV(2);
		param[3] += deltaV(3);
		param[4] += deltaV(4);
		param[5] += deltaV(5);

		if (deltaV.norm() < eps)
		{
			break;
		}
	}
}


template<class T>
void  __declspec(dllexport) Affine3dFittingQR(T param[9], T(&ortMat)[9], T(&tVec)[3]
	, T(*point0)[3], T(*point)[3], int point_num, double eps, int jt)
{
	VectorXd deltaV = VectorXd::Zero(9);
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

		param[0] += deltaV(0);
		param[1] += deltaV(1);
		param[2] += deltaV(2);
		param[3] += deltaV(3);
		param[4] += deltaV(4);
		param[5] += deltaV(5);
		param[6] += deltaV(6);
		param[7] += deltaV(7);
		param[8] += deltaV(8);

		if (deltaV.norm() < eps)
		{
			break;
		}
	}
}



template<class T>
void  __declspec(dllexport) AffineFittingQR(T param[5], T(&ortMat)[9]
	, T(*point0)[3], T(*point)[3], int point_num, double eps, int jt)
{
	VectorXd deltaV = VectorXd::Zero(5);
	param[0] = 0;
	param[1] = 0;
	param[2] = 1;
	param[3] = 0;
	param[4] = 1;

	for (int i = 0; i != jt; ++i)
	{
		GetDeltaVec4AffineQR(deltaV, ortMat, param, point0, point, point_num);

		param[0] += deltaV(0);
		param[1] += deltaV(1);
		param[2] += deltaV(2);
		param[3] += deltaV(3);
		param[4] += deltaV(4);

		if (deltaV.norm() < eps)
		{
			break;
		}
	}
}

template<class T>
void  __declspec(dllexport) NominalFitting(T nomiVal[2], T(*mat_set)[MAXMATSIZE], int deg
	, T(*point0)[3], T(*pointAff)[3], int point_num, double eps, int jt)
{
	VectorXd nomiDelta = VectorXd::Zero(2);
	T nomi_error = 0;
	for (int i = 0; i != jt; ++i)
	{

		nomi_error = GetDeltaVec4Nominal(nomiDelta, mat_set
			, nomiVal, deg, point0, pointAff, point_num);

		for (int j = 0; j != 2; ++j)
		{
			nomiVal[j] += nomiDelta(j);
		}

		//if (nomi_error < eps)
		//{
		//	break;
		//}
	}
}



inline void Analytic3dFitting(VectorXd &param, int param_num, int deg
	, double(*point0)[3], double(*pointAff)[3], int point_num, double eps, int jt)
{

	VectorXd deltaV = VectorXd::Zero(param_num);

	VectorXd hk = VectorXd::Zero(param_num);
	MatrixXd Hk = MatrixXd::Identity(param_num, param_num);
	VectorXd Xi = VectorXd::Zero(param_num);
	VectorXd delta_fi = VectorXd::Zero(param_num);

	VectorXd delta_Xi = VectorXd::Zero(param_num);
	VectorXd delta_gi = VectorXd::Zero(param_num);

	MatrixXd Uk = MatrixXd::Zero(param_num, param_num);

	double mat_error = 0;
	int sn_temp = 0;
	int temp;
	for (int i = 0; i != jt; ++i)
	{
		Hk = Hk + Uk;
		mat_error = GetDeltaVec3d4Analytic(deltaV, param_num, param
			, deg, point0, pointAff, point_num);

		hk = Hk*deltaV;
		sn_temp = 0;
		for (int j = 0; j <= deg; ++j)
		{
			temp = Combination(j + 2, 2) * 3;

			for (int k = 0; k != temp; ++k)
			{
				param(sn_temp + k) += hk(sn_temp + k);
			}

			sn_temp += temp;
		}

		if (mat_error < eps)
		{
			break;
		}

		delta_Xi = param - Xi;
		delta_gi = deltaV - delta_fi;

		BFGS(Uk, Hk, delta_gi, delta_Xi);

		Xi = param;
		delta_fi = deltaV;
	}
}


template<class T>
void  __declspec(dllexport) AnalyticFitting(T(*mat_set)[MAXMATSIZE], T nomiVal[3], int deg
	, T(*point0)[3], T(*pointAff)[3], int point_num, double eps, int jt)
{
	typedef T Td[MAXMATSIZE];
	Td *delta_set = new Td[MAXDEGREE];
	memset(delta_set, 0, sizeof(Td)*MAXDEGREE);
	T mat_error = 0;
	for (int i = 0; i != jt; ++i)
	{
		mat_error = GetDeltaVec4Analytic(delta_set, mat_set
			, nomiVal, deg, point0, pointAff, point_num);

		for (int j = 0; j <= deg; ++j)
		{
			for (int k = 0; k < 2 * (j + 1); ++k)
			{
				mat_set[j][k] += delta_set[j][k];
			}
		}

		if (mat_error < eps)
		{
			break;
		}
	}

	delete[] delta_set;
}

inline double CalcT(VectorXd &X0, VectorXd &hk, MatrixXd &matB, VectorXd &l)
{
	VectorXd bn = l - matB*X0;
	VectorXd fh0 = matB*hk;

	double a2 = (fh0.transpose()*fh0)(0);
	double ab = (fh0.transpose()*bn)(0);
	return ab / a2;
}

template<class T>
void  __declspec(dllexport) Quad3dFitting(VectorXd &param
	, T(*point0)[3], T(*pointAff)[3], int point_num, double eps, int jt)
{
	VectorXd deltaV = VectorXd::Zero(30);
	VectorXd hk = VectorXd::Zero(30);
	MatrixXd Hk = MatrixXd::Identity(30, 30);
	VectorXd Xi = VectorXd::Zero(30);
	VectorXd delta_fi = VectorXd::Zero(30);

	VectorXd delta_Xi = VectorXd::Zero(30);
	VectorXd delta_gi = VectorXd::Zero(30);

	MatrixXd Uk = MatrixXd::Zero(30, 30);
	T preNorm = INFINITYNUM;
	T rate = 1.0;
	for (int i = 0; i != jt; ++i)
	{
		Hk = Hk + Uk;
		MatrixXd matB = MatrixXd::Zero(3 * point_num, 30);
		VectorXd l = VectorXd::Zero(point_num * 3);
		GetDeltaVec4Quad3d(deltaV, matB, l, param, point0, pointAff, point_num);
		hk = Hk*deltaV;
		//rate = CalcT(param, hk, matB, l);
		param(0) += hk(0)*rate;
		param(1) += hk(1)*rate;
		param(2) += hk(2)*rate;
		param(3) += hk(3)*rate;
		param(4) += hk(4)*rate;
		param(5) += hk(5)*rate;
		param(6) += hk(6)*rate;
		param(7) += hk(7)*rate;
		param(8) += hk(8)*rate;
		param(9) += hk(9)*rate;
		param(10) += hk(10)*rate;
		param(11) += hk(11)*rate;
		param(12) += hk(12)*rate;
		param(13) += hk(13)*rate;
		param(14) += hk(14)*rate;
		param(15) += hk(15)*rate;
		param(16) += hk(16)*rate;
		param(17) += hk(17)*rate;
		param(18) += hk(18)*rate;
		param(19) += hk(19)*rate;
		param(20) += hk(20)*rate;
		param(21) += hk(21)*rate;
		param(22) += hk(22)*rate;
		param(23) += hk(23)*rate;
		param(24) += hk(24)*rate;
		param(25) += hk(25)*rate;
		param(26) += hk(26)*rate;
		param(27) += hk(27)*rate;
		param(28) += hk(28)*rate;
		param(29) += hk(29)*rate;
		T norm = deltaV.norm();

		if (norm < eps || norm>preNorm)
		{
			break;
		}
		preNorm = norm;

		delta_Xi = param - Xi;
		delta_gi = deltaV - delta_fi;

		BFGS(Uk, Hk, delta_gi, delta_Xi);

		Xi = param;
		delta_fi = deltaV;
	}
}


template<class T>
void __declspec(dllexport) GetAffine3dParamQR(T(&mat)[9], T(&t_a)[3]
	, T(&ortMat)[9], T(&tVec)[3], T(*ps0)[3], T(*ps)[3], int num)
{
	T param[9] = { 0 };
	Affine3dFittingQR(param, ortMat, tVec, ps0, ps, num, 0.01, 1);

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
	AffineFittingQR(param, ortMat, ps0, ps, num, 0.01, 1);

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
void __declspec(dllexport) GetAnalyTransParam(T(*mat_set)[MAXMATSIZE], T nomiVal[2], int deg
	, T(*ps0)[3], T(*ps)[3], int num)
{
	//ADMMFitting(mat_set, deg, ps0, ps, num, 0.01, 50);

	//for (int i = 2; i <= deg; ++i)
	{
		AnalyticFitting(mat_set, nomiVal, deg, ps0, ps, num, 0.00000001, 1);
		NominalFitting(nomiVal, mat_set, deg, ps0, ps, num, 0.01, 1);
	}
}

template<class T>
void __declspec(dllexport) GetAnalyTransf3dParam(T(*mat_set)[MAXMATSIZE]
	, int deg, T(*ps0)[3], T(*ps)[3], int num)
{
	int param_num = 0;
	for (int i = 0; i <= deg; ++i)
	{
		param_num += Combination(i + 2, 2) * 3;
	}
	VectorXd param = VectorXd::Zero(param_num);
	Analytic3dFitting(param, param_num, deg, ps0, ps, num, 0.00000001, 1);

	int sn_temp = 0;
	for (int j = 0; j <= deg; ++j)
	{
		int temp = Combination(j + 2, 2) * 3;

		for (int k = 0; k != temp; ++k)
		{
			mat_set[j][k] = param(sn_temp + k);
		}

		sn_temp += temp;
	}
}


inline double Analytic_ICP_3D(Vector3d *point_set4regist, int num4regist
	, Vector3d *aim_point_set, int aim_num, int icpIterCount, double eps)
{
	vector<int> fit_set;

	double dist;
	dist = 0;
	for (int i = 0; i != num4regist; ++i)
	{
		dist += sqrt(pow(aim_point_set[i](0) - point_set4regist[i](0), 2)
			+ pow(aim_point_set[i](1) - point_set4regist[i](1), 2)
			+ pow(aim_point_set[i](2) - point_set4regist[i](2), 2));

	}
	dist /= num4regist;

	printf("initial error: %.6f\n\n", dist);

	MatrixXd R = MatrixXd::Identity(3, 3);
	VectorXd t = VectorXd::Zero(3);
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

	int deg = 2;

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

		matrix3double3(d_set, point_set4regist, num4regist);

		double9 ort;
		double3 tV;

		ort[0] = R(0, 0);
		ort[1] = R(0, 1);
		ort[2] = R(0, 2);
		ort[3] = R(1, 0);
		ort[4] = R(1, 1);
		ort[5] = R(1, 2);
		ort[6] = R(2, 0);
		ort[7] = R(2, 1);
		ort[8] = R(2, 2);

		tV[0] = t(0);
		tV[1] = t(1);
		tV[2] = t(2);
		GetAffine3dParamQR(coef1, coef0, ort, tV, d_set, aim_d_set, num4regist);

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
			memcpy(d_set[i], d_t, sizeof(double3));

			//dist += get_euclid_dist<double, 3>(d_t, aim_d_set[i]);
			dist += (aim_point_set[i] - point_set4regist[i]).norm();

		}
		dist /= num4regist;

		if (i >= 1)
		{

			GetAnalyTransf3dParam(mat_set, deg, d_set, aim_d_set, num4regist);

			for (int i = 0; i != num4regist; ++i)
			{
				Analyt3dTrans(mat_set, deg, d_set[i], d_t);
				double2matrix(point_set4regist[i], d_t);

				//dist += get_euclid_dist<double, 2>(d_t, aim_d_set[i]);
				dist += (aim_point_set[i] - point_set4regist[i]).norm();
			}

			if (i % 2 == 0)
			{
				deg += 1;
			}

			dist /= num4regist;
			printf("-------------series error is %.15f-------------\n", (float)dist);
		}
		if (icpICount > 2)
		{
			icpICount /= 1.5;
		}

		if (dist > preDist && dist < eps)
		{
			result = dist;

			printf("-------------final error is %.15f-------------\n", (float)dist);
			break;
		}

		preDist = dist;

		finishStep = clock();
		durationStep = (double)(finishStep - start) / CLOCKS_PER_SEC;
		fprintf(fp, "%.15f,%f\n", (float)dist, durationStep);


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
inline double Analytic_ICP_2D(Vector2d *point_set4regist, int num4regist
	, Vector2d *aim_point_set, int aim_num, int icpIterCount, double eps)
{
	vector<int> fit_set;
	coefmat *mat_set = new coefmat[MAXDEGREE];
	MatrixXd R = MatrixXd::Identity(2, 2);
	VectorXd t = VectorXd::Zero(2);
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
	double3 nomiVal;
	memset(nomiVal, 0, sizeof(double3));

	FILE *fp = fopen("F:\\icp_error_time.csv", "w+");
	if (fp == NULL) {
		fprintf(stderr, "fopen() failed.\n");
		exit(EXIT_FAILURE);
	}

	clock_t start, finish, icpStart, icpFinish,finishStep;

	double duration, durationStep, icpDuration;

	start = clock();

	int deg = 2;

	for (int i = 0; i != MAXNEWICPITERCOUNT; ++i)
	{
		icpStart = clock();
		icp.setMaxIterations(icpICount);
		icp.fit(point_set4regist, num4regist, R, t, fit_set, -1);

		for (int j = 0; j != num4regist; ++j)
		{
			point_set4regist[j] = R*point_set4regist[j] + t;
		}

		for (int j = 0; j != fit_set.size(); ++j)
		{
			matrix2double3(aim_d_set[j], aim_point_set[fit_set[j]]);
		}

		if (0 == i)
		{
			rigidErr = 0;
			for (int j = 0; j != num4regist; ++j)
			{
				//matrix2double3(d_t, point_set4regist[j]);
				//rigidErr += get_euclid_dist<double, 2>(d_t, aim_d_set[j]);

				rigidErr += (aim_point_set[j] - point_set4regist[j]).norm();
			}

			rigidErr /= num4regist;

			printf("-------------rigid error is %f-------------\n", (float)rigidErr);

			icpFinish = clock();
			icpDuration = (double)(icpFinish - icpStart) / CLOCKS_PER_SEC;

			printf("icp running time: %.5f s\n", icpDuration);
		}



		matrix2double3(d_set, point_set4regist, num4regist);

		double9 ort;
		memset(ort, 0, sizeof(double9));
		ort[0] = R(0, 0);
		ort[1] = R(0, 1);
		ort[2] = t(0);
		ort[3] = R(1, 0);
		ort[4] = R(1, 1);
		ort[5] = t(1);
		ort[8] = 1;
		GetAffineTransformParamQR(matrix, ort, d_set, aim_d_set, num4regist);
		//GetPerspectiveTransform(matrix, d_set, aim_d_set, num4regist);
		//GetQuadraticTransformParam(matQ, matrix, d_set, aim_d_set, num4regist);
		//GetCubicTransformParam(matC, matQ, matrix, d_set, aim_d_set, num4regist);
		//memset(matQ, 0, sizeof(double9));

		dist = 0;
		memset(mat_set, 0, sizeof(coefmat)*MAXDEGREE);

		memset(nomiVal, 0, sizeof(double3));
		mat_set[0][0] = matrix[2];
		mat_set[0][1] = matrix[5];

		mat_set[1][0] = matrix[0];
		mat_set[1][1] = matrix[1];
		mat_set[1][2] = matrix[3];
		mat_set[1][3] = matrix[4];

		dist = 0;
		for (int i = 0; i != num4regist; ++i)
		{
			AnalyticTrans(mat_set, nomiVal, deg, d_set[i], d_t);
			double2matrix(point_set4regist[i], d_t);
			memcpy(d_set[i], d_t, sizeof(double3));

			dist += (aim_point_set[i] - point_set4regist[i]).norm();

		}
		dist /= num4regist;

		printf("-------------affine error is %.15f-------------\n", (float)dist);

		if (i >= 0)
		{

			GetAnalyTransParam(mat_set, nomiVal, deg, d_set, aim_d_set, num4regist);

			printf("%f\t%f\t\n%f\t%f\t%f\t%f\n\n"
				, mat_set[0][0], mat_set[0][1]
				, mat_set[1][0], mat_set[1][1]
				, mat_set[1][2], mat_set[1][3]);

			printf("nominal value: %f\t%f\n\n"
				, nomiVal[0], nomiVal[1]);


			for (int i = 0; i != num4regist; ++i)
			{
				AnalyticTrans(mat_set, nomiVal, deg, d_set[i], d_t);
				double2matrix(point_set4regist[i], d_t);

				//dist += get_euclid_dist<double, 2>(d_t, aim_d_set[i]);

				dist += (aim_point_set[i] - point_set4regist[i]).norm();
			}

			deg += 1;

			dist /= num4regist;
			printf("-------------series error is %.15f-------------\n", (float)dist);
		}
		if (icpICount > 2)
		{
			icpICount /= 1.5;
		}

		if (dist > preDist && dist < eps)
		{
			result = dist;

			printf("-------------final error is %.15f-------------\n", (float)dist);
			break;
		}

		preDist = dist;

		finishStep = clock();
		durationStep = (double)(finishStep - start) / CLOCKS_PER_SEC;
		fprintf(fp, "%.15f,%f\n", (float)dist, durationStep);
	}
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf("ftransform algorithms' response time is--------------%.5f seconds-------------", duration);

	fclose(fp);

	delete[]aim_d_set;
	delete[]mat_set;
	delete[] d_set;
	return result;
}


#endif