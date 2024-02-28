// dllmain.cpp : 定义 DLL 应用程序的入口点。
#include "stdafx.h"
#include "algorithm\Fitting.h"

BOOL APIENTRY DllMain(HMODULE hModule,
	DWORD  ul_reason_for_call,
	LPVOID lpReserved
	)
{
	switch (ul_reason_for_call)
	{
	case DLL_PROCESS_ATTACH:
	case DLL_THREAD_ATTACH:
	case DLL_THREAD_DETACH:
	case DLL_PROCESS_DETACH:
		break;
	}
	return TRUE;
}

inline double feature_nonrigid_transform(double(*point_set4regist)[2], int num4regist
	, double(*aim_point_set)[2], int aim_num, int iterCount)
{
	CPoint2 *ps4regist = new CPoint2[num4regist];
	CPoint2 *ps = new CPoint2[aim_num];

	for (int i = 0; i != num4regist; ++i)
	{
		ps4regist[i][0] = point_set4regist[i][0];
		ps4regist[i][1] = point_set4regist[i][1];
	}

	for (int i = 0; i != aim_num; ++i)
	{
		ps[i][0] = aim_point_set[i][0];
		ps[i][1] = aim_point_set[i][1];
	}

	double result = Analytic_ICP_2D(ps4regist, num4regist
		, ps, aim_num, iterCount, 0.01);

	for (int i = 0; i != num4regist; ++i)
	{
		point_set4regist[i][0] = ps4regist[i][0];
		point_set4regist[i][1] = ps4regist[i][1];
	}

	for (int i = 0; i != aim_num; ++i)
	{
		aim_point_set[i][0] = ps[i][0];
		aim_point_set[i][1] = ps[i][1];
	}

	delete[] ps4regist;
	delete[] ps;

	return result;
}

inline double feature_nonrigid_transf3d(double(*point_set4regist)[3], int num4regist
	, double(*aim_point_set)[3], int aim_num, int iterCount, double eps)
{
	CPoint3 *ps4regist = new CPoint3[num4regist];
	CPoint3 *ps = new CPoint3[aim_num];

	for (int i = 0; i != num4regist; ++i)
	{
		ps4regist[i][0] = point_set4regist[i][0];
		ps4regist[i][1] = point_set4regist[i][1];
		ps4regist[i][2] = point_set4regist[i][2];
	}

	for (int i = 0; i != aim_num; ++i)
	{
		ps[i][0] = aim_point_set[i][0];
		ps[i][1] = aim_point_set[i][1];
		ps[i][2] = aim_point_set[i][2];
	}

	double result = Analytic_ICP_3D(ps4regist, num4regist
		, ps, aim_num, iterCount, eps);

	for (int i = 0; i != num4regist; ++i)
	{
		point_set4regist[i][0] = ps4regist[i][0];
		point_set4regist[i][1] = ps4regist[i][1];
		point_set4regist[i][2] = ps4regist[i][2];
	}

	for (int i = 0; i != aim_num; ++i)
	{
		aim_point_set[i][0] = ps[i][0];
		aim_point_set[i][1] = ps[i][1];
		aim_point_set[i][2] = ps[i][2];
	}

	delete[] ps4regist;
	delete[] ps;

	return result;
}


DLL_ANALYTICICP double feature_transform(double(*point_set4regist)[2], int num4regist
	, double(*aim_point_set)[2], int aim_num, int iterCount)
{
	return feature_nonrigid_transform(point_set4regist, num4regist
		, aim_point_set, aim_num, iterCount);
}

DLL_ANALYTICICP double feature_transf3d(double(*point_set4regist)[3], int num4regist
	, double(*aim_point_set)[3], int aim_num, int iterCount, double eps)
{
	return feature_nonrigid_transf3d(point_set4regist, num4regist
		, aim_point_set, aim_num, iterCount, eps);
}

DLL_ANALYTICICP void randomly_analytic_map(double(*taget)[2], double(*moved)[2]
	, int num, int deg)
{
	AnalyticTransf(taget, moved, num, deg);
}

DLL_ANALYTICICP void randomly_3d_map(double(*taget)[3], double(*moved)[3]
	, int num)
{
	AnalyticT3d(taget, moved, num);
}