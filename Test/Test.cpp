// test.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include <windows.h>
#include <opencv2/opencv.hpp>
#include "time.h"
using namespace cv;
using namespace std;


int main()
{
	typedef bool(*ma_initialize)(TCHAR *dllPath);

	typedef double(*ma_ps_regist)(const char *fixed_path, const char *moving_path
		, unsigned char *bmp_qst, int r, int c);

	typedef void(*ma_analytic_map)(const char *fixed_path
		, unsigned char *bmp_qst, int r, int c);
	typedef void(*ma_release)();

	HINSTANCE hDll;
	printf("beginning----------------------------\n\n");
	hDll = LoadLibrary(L"SmoothAdjustment.dll");

	ma_initialize initialize = (ma_initialize)GetProcAddress(hDll, "carving_initialize");
	//2D Analytic-ICP
	//ma_ps_regist ps_regist = (ma_ps_regist)GetProcAddress(hDll, "points_regist");

	//3D second-order Analytic-ICP
	ma_ps_regist ps_regist = (ma_ps_regist)GetProcAddress(hDll, "p3d_regist");

	//2D analytic mapping
	//ma_analytic_map analytic_map = (ma_analytic_map)GetProcAddress(hDll, "analytic_points_map");

	//3D second-order analytic mapping
	ma_analytic_map analytic_map = (ma_analytic_map)GetProcAddress(hDll, "analytic_3d_map");
	ma_release release = (ma_release)GetProcAddress(hDll, "laser_release");

	int registType = 0;
	int addPerturb = 0;
	char buf[MAX_PATH] = { 0 };
	char movingPsPath[MAX_PATH] = { 0 };
	char fixedPsPath[MAX_PATH] = { 0 };
	char outputPerturbImgPath[MAX_PATH] = { 0 };
	char outputEffectImgPath[MAX_PATH] = { 0 };

	int len = GetPrivateProfileStringA("Param",
		"RegistType",
		"0",
		buf,
		MAX_PATH,
		".\\experiment.ini");

	registType = atoi(buf);

	len = GetPrivateProfileStringA("Param",
		"AddPerturb",
		"0",
		buf,
		MAX_PATH,
		".\\experiment.ini");

	addPerturb = atoi(buf);

	len = GetPrivateProfileStringA("Param",
		"MovingPsPath",
		"",
		movingPsPath,
		MAX_PATH,
		".\\experiment.ini");

	len = GetPrivateProfileStringA("Param",
		"FixedPsPath",
		"",
		fixedPsPath,
		MAX_PATH,
		".\\experiment.ini");

	len = GetPrivateProfileStringA("Param",
		"OutputPerturbImgPath",
		"",
		outputPerturbImgPath,
		MAX_PATH,
		".\\experiment.ini");

	len = GetPrivateProfileStringA("Param",
		"OutputEffectImgPath",
		"",
		outputEffectImgPath,
		MAX_PATH,
		".\\experiment.ini");

	initialize(_T("Analytic_ICP.dll"));

	printf("Non-rigid registration with analytic mapping begin------------\n");

	Mat mat;
	mat.create(2048, 2048, CV_8UC3);

	if (1 == registType)
	{
		if (addPerturb)
		{
			analytic_map(movingPsPath
				, mat.data, mat.rows, mat.cols);
			imwrite(outputPerturbImgPath, mat);
		}

		//3D
		double result = ps_regist(fixedPsPath
			, movingPsPath
			, mat.data, mat.rows, mat.cols);
		imwrite(outputEffectImgPath, mat);
	}
	else
	{
		//2D
		if (addPerturb)
		{
			analytic_map(movingPsPath
				, mat.data, mat.rows, mat.cols);
			imwrite(outputPerturbImgPath, mat);
		}

		//2D
		double result = ps_regist(fixedPsPath
			, movingPsPath
			, mat.data, mat.rows, mat.cols);
		imwrite(outputEffectImgPath, mat);
	}

	release();
	system("pause");
	FreeLibrary(hDll);
	return 0;
}
