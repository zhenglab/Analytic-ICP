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

	initialize(_T("Analytic_ICP.dll"));

	printf("Non-rigid registration with analytic mapping begin------------\n");

	Mat mat;
	mat.create(2048, 2048, CV_8UC3);

	//3D
	analytic_map("F:\\bunny test\\bunny_3500.csv"
		, mat.data, mat.rows, mat.cols);
	imwrite("F:\\analymap.bmp", mat);

	//3D
	double result = ps_regist("F:\\fixed3d.csv"
		, "F:\\bunny test\\bunny_3500.csv"
		, mat.data, mat.rows, mat.cols);
	imwrite("F:\\ps_regist_test.bmp", mat);


	////2D
	//analytic_map("F:\\fixtures\\fish.csv"
	//	, mat.data, mat.rows, mat.cols);
	//imwrite("F:\\analymap.bmp", mat);

	////2D
	//double result = ps_regist("C:\\Users\\Administrator\\Desktop\\2ddeformandresidual\\0.06\\fixed_2.csv"
	//	, "F:\\fixtures\\fish.csv"
	//	, mat.data, mat.rows, mat.cols);
	//imwrite("F:\\ps_regist_test.bmp", mat);

	release();
	system("pause");
	FreeLibrary(hDll);
	return 0;
}