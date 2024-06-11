/*!
*      \brief the entry for interface
*	   \author Wei Feng
*      \date 09/16/2020
*/
#include "stdafx.h"
#include "algorithms\Algo.h"
#include "fileM\FileProcess.h"


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

DLL_PUREENTITY bool carving_initialize(TCHAR *dllPath)
{
	hDll = LoadLibrary(dllPath);

	ftransform = (Ftransform)GetProcAddress(hDll, "feature_transform");
	ftransf3d = (Ftransf3d)GetProcAddress(hDll, "feature_transf3d");
	analyMap = (RandomAnalyM)GetProcAddress(hDll, "randomly_analytic_map");
	analyMap3d = (Random3dM)GetProcAddress(hDll, "randomly_3d_map");
	return true;
}


/*!
free the memory
*/
DLL_PUREENTITY void laser_release()
{
	FreeLibrary(hDll);
}


DLL_PUREENTITY double p3d_regist(const char *fixed_path, const char *moving_path
	, const char *moving_output)
{
	double3 *fixed = new double3[MAX_VERTEX_NUM];
	int *fixed_idx = new int[MAX_VERTEX_NUM];
	int f_num = 0;

	double3 *moving_prev = new double3[MAX_VERTEX_NUM];
	double3 *moving = new double3[MAX_VERTEX_NUM];
	int *mov_idx = new int[MAX_VERTEX_NUM];
	int m_num = 0;

	data_from_path(fixed, f_num, fixed_path, 1.0);
	data_from_path(moving, m_num, moving_path, 1.0);

	memcpy(moving_prev, moving, sizeof(double3)*m_num);

	double result;
	result = ftransf3d(moving, m_num, fixed, f_num, 3);

	write2csv(moving, m_num, "F:\\moved.csv");

	delete[]fixed;
	delete[]moving;
	delete[]moving_prev;
	return result;
}



DLL_PUREENTITY double points_regist(const char *fixed_path, const char *moving_path
	, unsigned char *bmp_qst, int r, int c)
{
	double2 *fixed = new double2[MAX_VERTEX_NUM];
	int f_num = 0;

	double2 *moving_prev = new double2[MAX_VERTEX_NUM];
	double2 *moving = new double2[MAX_VERTEX_NUM];
	int m_num = 0;

	data_from_path(fixed, f_num, fixed_path, 1.0);
	data_from_path(moving, m_num, moving_path, 1.0);

	memcpy(moving_prev, moving, sizeof(double2)*m_num);
	double result;
	result = ftransform(moving, m_num, fixed, f_num, 2);

	double2 *merging = new double2[f_num + m_num + m_num];

	memcpy(merging, fixed, sizeof(double2)*f_num);
	memcpy(merging + f_num, moving, sizeof(double2)*m_num);
	memcpy(merging + f_num + m_num, moving_prev, sizeof(double2)*m_num);

	for (int i = 0; i != f_num + m_num + m_num; ++i)
	{
		merging[i][0] *= 300;
		merging[i][1] *= 300;
	}

	Rect boundingbox;
	getboundingbox(boundingbox, merging, f_num + m_num + m_num, 0);

	for (int i = 0; i != f_num + m_num + m_num; ++i)
	{
		merging[i][0] -= boundingbox.x - 150;
		merging[i][1] -= boundingbox.y - 150;
	}

	for (int i = f_num + m_num; i != f_num + m_num + m_num; ++i)
	{
		merging[i][0] += 900;
	}

	Mat mat;
	mat.create(r, c, CV_8UC3);
	memset(mat.data, 0, sizeof(uchar) * 3 * r*c);
	for (int i = 0; i != r; ++i)
	{
		for (int j = 0; j != c; ++j)
		{
			cv::Vec3b temp_channel_array;
			temp_channel_array[0] = 255;
			temp_channel_array[1] = 255;
			temp_channel_array[2] = 255;
			mat.at<cv::Vec3b>(i, j) = temp_channel_array;
		}
	}
	Data2bmp(mat, merging, f_num, 0, 0, 255);

	Data2bmp(mat, merging + f_num, m_num, 0, 255, 0);

	Data2bmp(mat, merging + f_num + m_num, m_num, 255, 0, 0);

	memcpy(bmp_qst, mat.data, sizeof(uchar) * 3 * r * c);

	return result;
}


DLL_PUREENTITY void analytic_points_map(const char *fixed_path
	, unsigned char *bmp_qst, int r, int c)
{
	double2 *moving = new double2[MAX_VERTEX_NUM];
	int f_num = 0;

	double2 *fixed = new double2[MAX_VERTEX_NUM];
	double2 *fixed_c = new double2[MAX_VERTEX_NUM];

	data_from_path(moving, f_num, fixed_path, 1.0);

	int deg = 8;
	Mat mat;
	mat.create(r, c, CV_8UC3);

	analyMap(moving, fixed, f_num, deg);

	write2csv(fixed, f_num, "F:\\fixed.csv");
	for (int i = 0; i != f_num; ++i)
	{
		fixed_c[i][0] = 500 * fixed[i][0] + 1024;
		fixed_c[i][1] = 500 * fixed[i][1] + 800;
	}

	//memset(mat.data, 0, sizeof(uchar) * 3 * r*c);
	for (int i = 0; i != 3 * r*c; ++i)
	{
		mat.data[i] = 255;
	}
	Data2bmp(mat, fixed_c, f_num, 0, 0, 0);

	memcpy(bmp_qst, mat.data, sizeof(uchar) * 3 * r * c);

	delete[]moving;
	delete[]fixed;
	delete[]fixed_c;
}

DLL_PUREENTITY void analytic_3d_map(const char *fixed_path
	, unsigned char *bmp_qst, int r, int c)
{
	double3 *moving = new double3[MAX_VERTEX_NUM];

	memset(moving, 0, sizeof(double3)*MAX_VERTEX_NUM);
	int ids[MAX_VERTEX_NUM] = { 0 };
	int f_num = 0;

	double3 *fixed = new double3[MAX_VERTEX_NUM];
	memset(fixed, 0, sizeof(double3)*MAX_VERTEX_NUM);
	double3 *fixed_c = new double3[MAX_VERTEX_NUM];
	memset(fixed_c, 0, sizeof(double3)*MAX_VERTEX_NUM);
	data_from_path(moving, f_num, fixed_path, 1.0);
	write2csv(moving, f_num, "F:\\moving3d.csv");
	int deg = 2;
	Mat mat;
	mat.create(r, c, CV_8UC3);

	analyMap3d(moving, fixed, f_num, 2);

	write2m(fixed, ids, f_num, "F:\\fixed3d.m");
	write2csv(fixed, f_num, "F:\\fixed3d.csv");
	for (int i = 0; i != f_num; ++i)
	{
		fixed_c[i][0] = 500 * fixed[i][0] + 1024;
		fixed_c[i][1] = 500 * fixed[i][1] + 300;
		fixed_c[i][2] = 500 * fixed[i][2] + 100;
	}

	memset(mat.data, 0, sizeof(uchar) * 3 * r*c);
	Data2bmp(mat, fixed_c, f_num, 255, 255, 255);

	memcpy(bmp_qst, mat.data, sizeof(uchar) * 3 * r * c);
}