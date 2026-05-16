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

DLL_PUREENTITY bool aicp_initialize(TCHAR* dllPath)
{
	hDll = LoadLibrary(dllPath);
	if (!hDll)
	{
		return false;
	}

	ftransform = (Ftransform)GetProcAddress(hDll, "feature_transf2d");
	ftransf3d = (Ftransf3d)GetProcAddress(hDll, "feature_transf3d");
	analyMap = (RandomAnalyM)GetProcAddress(hDll, "randomly_2d_map");
	analyMap3d = (Random3dM)GetProcAddress(hDll, "randomly_3d_map");

	if (!ftransform || !ftransf3d || !analyMap || !analyMap3d)
	{
		FreeLibrary(hDll);
		hDll = NULL;
		return false;
	}

	return true;
}


DLL_PUREENTITY void aicp_release()
{
	if (hDll)
	{
		FreeLibrary(hDll);
		hDll = NULL;
	}

	ftransform = NULL;
	ftransf3d = NULL;
	analyMap = NULL;
	analyMap3d = NULL;
}

DLL_PUREENTITY double p3d_regist(
	const char* fixed_path,
	const char* moving_path,
	const char* moved_output)
{
	if (!fixed_path || !moving_path || !moved_output)
	{
		return -1.0;
	}

	if (!ftransf3d)
	{
		return -1.0;
	}

	double3* fixed = new double3[MAX_VERTEX_NUM];
	int f_num = 0;

	double3* moving = new double3[MAX_VERTEX_NUM];
	int m_num = 0;

	data_from_path(fixed, f_num, fixed_path, 1.0);
	data_from_path(moving, m_num, moving_path, 1.0);

	if (f_num <= 0 || m_num <= 0)
	{
		delete[] fixed;
		delete[] moving;
		return -1.0;
	}

	double result = ftransf3d(moving, m_num, fixed, f_num, 3);

	write2csv(moving, m_num, moved_output);

	delete[] fixed;
	delete[] moving;

	return result;
}



DLL_PUREENTITY double p2d_regist(const char *fixed_path, const char *moving_path
	, const char *moved_output)
{

	if (!fixed_path || !moving_path || !moved_output)
	{
		return -1.0;
	}

	if (!ftransf3d)
	{
		return -1.0;
	}

	double2 *fixed = new double2[MAX_VERTEX_NUM];
	int f_num = 0;

	double2 *moving = new double2[MAX_VERTEX_NUM];
	int m_num = 0;

	data_from_path(fixed, f_num, fixed_path, 1.0);
	data_from_path(moving, m_num, moving_path, 1.0);


	if (f_num <= 0 || m_num <= 0)
	{
		delete[] fixed;
		delete[] moving;
		return -1.0;
	}

	double result;
	result = ftransform(moving, m_num, fixed, f_num, 2);
	write2csv(moving, m_num, moved_output);

	delete[] fixed;
	delete[] moving;

	return result;
}


DLL_PUREENTITY void analytic_2d_map(const char *fixed_path
	, const char *deformed_output, int deg)
{
	double2 *moving = new double2[MAX_VERTEX_NUM];
	int f_num = 0;

	double2 *fixed = new double2[MAX_VERTEX_NUM];

	data_from_path(moving, f_num, fixed_path, 1.0);

	analyMap(moving, fixed, f_num, deg);
	write2csv(fixed, f_num, deformed_output);
	

	delete[]moving;
	delete[]fixed;
}

DLL_PUREENTITY void analytic_3d_map(const char *fixed_path
	, const char *deformed_output, int deg)
{
	double3 *moving = new double3[MAX_VERTEX_NUM];

	memset(moving, 0, sizeof(double3)*MAX_VERTEX_NUM);
	int f_num = 0;

	double3 *fixed = new double3[MAX_VERTEX_NUM];
	memset(fixed, 0, sizeof(double3)*MAX_VERTEX_NUM);
	data_from_path(moving, f_num, fixed_path, 1.0);

	analyMap3d(moving, fixed, f_num, deg);
	write2csv(fixed, f_num, deformed_output);

	delete[]moving;
	delete[]fixed;
}