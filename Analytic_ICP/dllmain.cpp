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

inline double feature_nonrigid_regist(double(*point_set4regist)[2], int num4regist
	, double(*aim_point_set)[2], int aim_num, int iterCount)
{
	double expect4regist[2] = { 0 };
	double vari4regist = 0;
	double expect4aim[2] = { 0 };
	double vari4aim = 0;
	ExpectandVari<double, 2>(expect4regist, vari4regist, point_set4regist, num4regist);
	ExpectandVari<double, 2>(expect4aim, vari4aim, aim_point_set, aim_num);
	double vari = max(vari4regist, vari4aim);
	Vector2d *ps4regist = new Vector2d[num4regist];
	Vector2d *ps = new Vector2d[aim_num];
	for (int i = 0; i != num4regist; ++i)
	{
		ps4regist[i](0) = (point_set4regist[i][0] - expect4regist[0]) / vari;
		ps4regist[i](1) = (point_set4regist[i][1] - expect4regist[1]) / vari;
	}

	for (int i = 0; i != aim_num; ++i)
	{
		ps[i](0) = (aim_point_set[i][0] - expect4aim[0]) / vari;
		ps[i](1) = (aim_point_set[i][1] - expect4aim[1]) / vari;
	}

	double result = Analytic_ICP_2D(ps4regist, num4regist
		, ps, aim_num, iterCount, 0.0001);

	for (int i = 0; i != num4regist; ++i)
	{
		point_set4regist[i][0] = ps4regist[i][0] * vari /*+ expect4regist[0]*/;
		point_set4regist[i][1] = ps4regist[i][1] * vari /*+ expect4regist[1]*/;
	}

	for (int i = 0; i != aim_num; ++i)
	{
		aim_point_set[i][0] = ps[i][0] * vari /*+ expect4aim[0]*/;
		aim_point_set[i][1] = ps[i][1] * vari/* + expect4aim[1]*/;
	}

	delete[] ps4regist;
	delete[] ps;

	return result;
}

inline double feature_nonrigid_regist3d(double(*point_set4regist)[3], int num4regist
	, double(*aim_point_set)[3], int aim_num, int iterCount)
{
	double expect4regist[3] = { 0 };
	double vari4regist = 0;
	double expect4aim[3] = { 0 };
	double vari4aim = 0;
	ExpectandVari<double, 3>(expect4regist, vari4regist, point_set4regist, num4regist);
	ExpectandVari<double, 3>(expect4aim, vari4aim, aim_point_set, aim_num);
	double vari = max(vari4regist, vari4aim);
	Vector3d *ps4regist = new Vector3d[num4regist];
	Vector3d *ps = new Vector3d[aim_num];
	for (int i = 0; i != num4regist; ++i)
	{
		ps4regist[i](0) = (point_set4regist[i][0] - expect4regist[0]) / vari;
		ps4regist[i](1) = (point_set4regist[i][1] - expect4regist[1]) / vari;
		ps4regist[i](2) = (point_set4regist[i][2] - expect4regist[2]) / vari;
	}

	for (int i = 0; i != aim_num; ++i)
	{
		ps[i](0) = (aim_point_set[i][0] - expect4aim[0]) / vari;
		ps[i](1) = (aim_point_set[i][1] - expect4aim[1]) / vari;
		ps[i](2) = (aim_point_set[i][2] - expect4aim[2]) / vari;
	}

	double result = Analytic_ICP_3D(ps4regist, num4regist
		, ps, aim_num, iterCount, 0.00000001);

	for (int i = 0; i != num4regist; ++i)
	{
		point_set4regist[i][0] = ps4regist[i][0] * vari + expect4aim[0];
		point_set4regist[i][1] = ps4regist[i][1] * vari + expect4aim[1];
		point_set4regist[i][2] = ps4regist[i][2] * vari + expect4aim[2];
	}

	//for (int i = 0; i != aim_num; ++i)
	//{
	//	aim_point_set[i][0] = ps[i][0] * vari - expect4aim[0];
	//	aim_point_set[i][1] = ps[i][1] * vari - expect4aim[1];
	//	aim_point_set[i][2] = ps[i][2] * vari - expect4aim[2];
	//}

	delete[] ps4regist;
	delete[] ps;

	return result;
}


DLL_ANALYTICICP double feature_transform(double(*point_set4regist)[2], int num4regist
	, double(*aim_point_set)[2], int aim_num, int iterCount)
{
	return feature_nonrigid_regist(point_set4regist, num4regist
		, aim_point_set, aim_num, iterCount);
}

DLL_ANALYTICICP double feature_transf3d(double(*point_set4regist)[3], int num4regist
	, double(*aim_point_set)[3], int aim_num, int iterCount)
{
	return feature_nonrigid_regist3d(point_set4regist, num4regist
		, aim_point_set, aim_num, iterCount);
}

DLL_ANALYTICICP void randomly_analytic_map(double(*taget)[2], double(*moved)[2]
	, int num, int deg)
{
	AnalyticTransf(taget, moved, num, deg);
}

DLL_ANALYTICICP void randomly_3d_map(double(*taget)[3], double(*moved)[3]
	, int num, int deg)
{
	AnalyticT3d(taget, moved, num,deg);
}