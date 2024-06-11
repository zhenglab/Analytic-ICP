// SmoothAdjustment.cpp : 定义 DLL 应用程序的导出函数。
//

#include "stdafx.h"
#include "SmoothAdjustment.h"


// 这是导出变量的一个示例
SMOOTHADJUSTMENT_API int nSmoothAdjustment=0;

// 这是导出函数的一个示例。
SMOOTHADJUSTMENT_API int fnSmoothAdjustment(void)
{
	return 42;
}

// 这是已导出类的构造函数。
// 有关类定义的信息，请参阅 SmoothAdjustment.h
CSmoothAdjustment::CSmoothAdjustment()
{
	return;
}
