#pragma once
/*!
*      \file FileProcess.h
*      \brief read and write files
*	   \author Wei Feng
*      \date 09/29/2018
*/

#include <time.h>
#include <direct.h>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <list>
#include "strutil.h"
#include "../algorithms/Algo.h"

#define MAX_LINE 2048

/*!
non-rigid registration algorithm interface
*/
typedef double(*Ftransform)(double(*point_set4regist)[2], int num4regist
	, double(*aim_point_set)[2], int aim_num, int iterCount);

/*!
3d rank 2 non-rigid registration algorithm interface
*/
typedef double(*Ftransf3d)(double(*point_set4regist)[3], int num4regist
	, double(*aim_point_set)[3], int aim_num, int iterCount);

/*!
a randomly generated analytic map interface
*/
typedef void(*RandomAnalyM)(double(*taget)[2], double(*moved)[2]
	, int num, int deg);

/*!
a randomly generated 3d rank 2 analytic map interface
*/
typedef void(*Random3dM)(double(*taget)[3], double(*moved)[3]
	, int num, int deg);

Ftransform ftransform;
Ftransf3d ftransf3d;
RandomAnalyM analyMap;
Random3dM analyMap3d;
HINSTANCE hDll;

/*! \brief Mess4PlugIn struct
*   parameters for plug in
*/
struct Param4PlugIn
{
	char name[MAX_PATH];
	char path[MAX_PATH];
};


void data3_from_path(double(*ds)[3], int *ids, int &num, const string& path)
{
	std::fstream is(path, std::fstream::in);

	if (is.fail())
	{
		fprintf(stderr, "Error in opening file %s\n", path);
		return;
	}

	char buffer[MAX_LINE];
	num = 0;

	while (is.getline(buffer, MAX_LINE))
	{

		std::string line(buffer);
		line = strutil::trim(line);

		strutil::Tokenizer stokenizer(line, " \r\n");

		stokenizer.nextToken();
		std::string token = stokenizer.getToken();

		if (token == "Vertex")
		{
			stokenizer.nextToken();
			token = stokenizer.getToken();
			ids[num] = strutil::parseString<int>(token);

			double3 p;
			for (int j = 0; j < 3; j++)
			{
				stokenizer.nextToken();
				token = stokenizer.getToken();
				p[j] = strutil::parseString<float>(token);
			}

			std::memcpy(ds[num], p, sizeof(double3));

			if (!stokenizer.nextToken("\t\r\n"))
			{
				num++;
				continue;
			}

			token = stokenizer.getToken();

			int sp = (int)token.find("{");
			int ep = (int)token.find("}");

			if (sp >= 0 && ep >= 0)
			{
				//v->string() = token.substr(sp + 1, ep - sp - 1);
			}
			num++;
			continue;
		}
	}
}



template<int N>
void _declspec(dllexport) data_from_path(double (*ds)[N], int &num, const string& path, double zoom)
{
	num = 0;
	ifstream file(path);
	if (!file.is_open()) 
	{
		stringstream msg;
		msg << "Unable to open file for reading: " << path;
		throw runtime_error(msg.str());
	}
	string line;
	int i;
	while (getline(file, line)) 
	{
		i = 0;
		stringstream ss(line);
		double c;
		while (ss >> c) 
		{
			ds[num][i++] = c*zoom;
			// TODO support other delimiters than commas
			if (ss.peek() == ',') 
			{
				ss.ignore();
			}
		}
		num++;
	}
}


inline void write2csv(double(*taget)[2], int num, string path)
{
	FILE *fp = fopen(path.c_str(), "w+");
	if (fp == NULL) {
		fprintf(stderr, "fopen() failed.\n");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i != num; ++i)
	{
		if (i>0)
			fprintf(fp, "\n");
		fprintf(fp, "%f,%f", taget[i][0], taget[i][1]);
	}
	fclose(fp);
}

inline void write2csv(double(*taget)[3], int num, string path)
{
	FILE *fp = fopen(path.c_str(), "w+");
	if (fp == NULL) {
		fprintf(stderr, "fopen() failed.\n");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i != num; ++i)
	{
		fprintf(fp, "%f,%f,%f\n", taget[i][0], taget[i][1], taget[i][2]);
	}
	fclose(fp);
}

inline void write2m(double(*taget)[3], int *ids, int num, string path)
{
	FILE *fp = fopen(path.c_str(), "w+");
	if (fp == NULL) {
		fprintf(stderr, "fopen() failed.\n");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i != num; ++i)
	{
		fprintf(fp, "Vertex %d %f %f %f\n"
			, ids[i], (float)taget[i][0], (float)taget[i][1], (float)taget[i][2]);
	}
	fclose(fp);
}

inline void Data2bmp(Mat &op, double2 *ps, int num
	, BYTE r, BYTE g, BYTE b)
{
	for (int i = 0; i != num; ++i)
	{
		circle(op, cv::Point(ps[i][0], ps[i][1]), 5, cv::Scalar(r, g, b), 2);
	}
}

inline void Data2bmp(Mat &op, double3 *ps, int num
	, BYTE r, BYTE g, BYTE b)
{
	for (int i = 0; i != num; ++i)
	{
		circle(op, cv::Point(ps[i][0], ps[i][1]), 5, cv::Scalar(r, g, b), 2);
	}
}
