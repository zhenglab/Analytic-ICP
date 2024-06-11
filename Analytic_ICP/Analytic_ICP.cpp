// Analytic_ICP.cpp : 定义 DLL 应用程序的导出函数。
//

#include "stdafx.h"
#include "Analytic_ICP.h"


// 这是导出变量的一个示例
ANALYTIC_ICP_API int nAnalytic_ICP=0;

// 这是导出函数的一个示例。
ANALYTIC_ICP_API int fnAnalytic_ICP(void)
{
	return 42;
}

// 这是已导出类的构造函数。
// 有关类定义的信息，请参阅 Analytic_ICP.h
CAnalytic_ICP::CAnalytic_ICP()
{
	return;
}
