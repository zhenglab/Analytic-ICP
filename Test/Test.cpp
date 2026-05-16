// test.cpp : Console frontend for Analytic-ICP / Structured Analytic Mappings
//

#include "stdafx.h"
#include <iostream>
#include <windows.h>
#include <tchar.h>
#include <cstdlib>
#include <cstring>

using namespace std;

static int ReadIniInt(
	const char* section,
	const char* key,
	int defaultValue,
	const char* iniPath)
{
	char buf[MAX_PATH] = { 0 };
	char defaultBuf[64] = { 0 };

	sprintf_s(defaultBuf, "%d", defaultValue);

	GetPrivateProfileStringA(
		section,
		key,
		defaultBuf,
		buf,
		MAX_PATH,
		iniPath);

	return atoi(buf);
}

static bool ReadIniString(
	const char* section,
	const char* key,
	const char* defaultValue,
	char* output,
	DWORD outputSize,
	const char* iniPath)
{
	DWORD len = GetPrivateProfileStringA(
		section,
		key,
		defaultValue,
		output,
		outputSize,
		iniPath);

	return len > 0;
}

static void PrintLoadLibraryError(const char* dllName)
{
	DWORD err = GetLastError();

	char currentDir[MAX_PATH] = { 0 };
	GetCurrentDirectoryA(MAX_PATH, currentDir);

	cerr << "Error: failed to load " << dllName << "." << endl;
	cerr << "Current working directory: " << currentDir << endl;
	cerr << "GetLastError = " << err << endl;

	LPVOID msgBuf = nullptr;
	FormatMessageA(
		FORMAT_MESSAGE_ALLOCATE_BUFFER |
		FORMAT_MESSAGE_FROM_SYSTEM |
		FORMAT_MESSAGE_IGNORE_INSERTS,
		NULL,
		err,
		MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
		(LPSTR)&msgBuf,
		0,
		NULL);

	if (msgBuf)
	{
		cerr << "System message: " << (char*)msgBuf << endl;
		LocalFree(msgBuf);
	}
}

int main()
{
	typedef bool(*ma_initialize)(TCHAR* dllPath);

	typedef double(*ma_ps_regist)(
		const char* fixed_path,
		const char* moving_path,
		const char* moved_output);

	typedef void(*ma_analytic_map)(
		const char* input_path,
		const char* deformed_output,
		int deg);

	typedef void(*ma_release)();

	const char* experimentIni = ".\\experiment.ini";

	printf("Structured Analytic Mappings / Analytic-ICP console frontend begins ----------------------------\n\n");

	// ---------------------------------------------------------------------
	// Load SmoothAdjustment.dll
	// ---------------------------------------------------------------------
	HINSTANCE hDll = LoadLibrary(L"SmoothAdjustment.dll");
	if (!hDll)
	{
		PrintLoadLibraryError("SmoothAdjustment.dll");
		cerr << "Please make sure SmoothAdjustment.dll and its dependent DLLs are in the runtime directory." << endl;
		system("pause");
		return -1;
	}

	ma_initialize initialize =
		(ma_initialize)GetProcAddress(hDll, "aicp_initialize");

	ma_release release =
		(ma_release)GetProcAddress(hDll, "aicp_release");

	ma_ps_regist regist2d =
		(ma_ps_regist)GetProcAddress(hDll, "p2d_regist");

	ma_ps_regist regist3d =
		(ma_ps_regist)GetProcAddress(hDll, "p3d_regist");

	ma_analytic_map analyticMap2d =
		(ma_analytic_map)GetProcAddress(hDll, "analytic_2d_map");

	ma_analytic_map analyticMap3d =
		(ma_analytic_map)GetProcAddress(hDll, "analytic_3d_map");

	if (!initialize || !release || !regist2d || !regist3d || !analyticMap2d || !analyticMap3d)
	{
		cerr << "Error: failed to load required functions from SmoothAdjustment.dll." << endl;
		FreeLibrary(hDll);
		system("pause");
		return -1;
	}

	// ---------------------------------------------------------------------
	// Read experiment.ini
	// ---------------------------------------------------------------------
	int registType = ReadIniInt("Param", "RegistType", 1, experimentIni);
	int addPerturb = ReadIniInt("Param", "AddPerturb", 0, experimentIni);
	int perturbDeg2D = ReadIniInt("Param", "PerturbDeg2D", 8, experimentIni);
	int perturbDeg3D = ReadIniInt("Param", "PerturbDeg3D", 2, experimentIni);

	char movingPsPath[MAX_PATH] = { 0 };
	char fixedPsPath[MAX_PATH] = { 0 };
	char perturbedOutputPath[MAX_PATH] = { 0 };
	char movedOutputPath[MAX_PATH] = { 0 };

	ReadIniString("Param", "MovingPsPath", "", movingPsPath, MAX_PATH, experimentIni);
	ReadIniString("Param", "FixedPsPath", "", fixedPsPath, MAX_PATH, experimentIni);
	ReadIniString("Param", "PerturbedOutputPath", "results\\perturbed.csv", perturbedOutputPath, MAX_PATH, experimentIni);
	ReadIniString("Param", "MovedOutputPath", "results\\moved.csv", movedOutputPath, MAX_PATH, experimentIni);

	if (strlen(movingPsPath) == 0 || strlen(fixedPsPath) == 0)
	{
		cerr << "Error: MovingPsPath or FixedPsPath is empty in experiment.ini." << endl;
		FreeLibrary(hDll);
		system("pause");
		return -1;
	}

	// ---------------------------------------------------------------------
	// Initialize Analytic_ICP.dll through SmoothAdjustment.dll
	// ---------------------------------------------------------------------
	if (!initialize(_T("Analytic_ICP.dll")))
	{
		cerr << "Error: failed to initialize Analytic_ICP.dll." << endl;
		FreeLibrary(hDll);
		system("pause");
		return -1;
	}

	printf("Analytic-ICP registration begins ------------------------------------------\n");
	printf("Fixed point set : %s\n", fixedPsPath);
	printf("Moving point set: %s\n", movingPsPath);
	printf("Output point set: %s\n", movedOutputPath);

	double result = -1.0;

	if (registType == 1)
	{
		// 3D
		if (addPerturb)
		{
			printf("Generating 3D analytic perturbation: degree = %d\n", perturbDeg3D);
			analyticMap3d(movingPsPath, perturbedOutputPath, perturbDeg3D);
			printf("Perturbed point set written to: %s\n", perturbedOutputPath);
		}

		result = regist3d(fixedPsPath, movingPsPath, movedOutputPath);
	}
	else
	{
		// 2D
		if (addPerturb)
		{
			printf("Generating 2D analytic perturbation: degree = %d\n", perturbDeg2D);
			analyticMap2d(movingPsPath, perturbedOutputPath, perturbDeg2D);
			printf("Perturbed point set written to: %s\n", perturbedOutputPath);
		}

		result = regist2d(fixedPsPath, movingPsPath, movedOutputPath);
	}

	printf("\nRegistration finished.\n");
	printf("Final error: %.10f\n", result);
	printf("Registered moving point set written to: %s\n", movedOutputPath);

	release();
	FreeLibrary(hDll);

	printf("Structured Analytic Mappings / Analytic-ICP console frontend finishes --------------------------\n");

	system("pause");
	return 0;
}