#pragma once
/*!
*      \file Algo.h
*      \brief algorithm
*      \author Wei Feng
*      \date 09/16/2020
*/
#ifndef MATHEMATICAL
#define MATHEMATICAL
#include <opencv2/opencv.hpp>
#include<math.h>
#include <vector>
#include <Windows.h>

#define OK 0
#define BREAK 1
#define COLORDIFF 2
#define SKEWSIZE 4
#define MULTICARVE 8
#define CHARNUMERROR 16
#define NONECODE 32
#define MISMATCH 64
#define SNCODENUMERROR 128
#define BITMAPERROR 256
#define INFINITYNUM 16777216
#define MAX_VERTEX_NUM 4096
#define BLACK 0
#define WHITE 255
#define GRAY 128
#define MAXITERCOUNT 10
#define PI 3.14159265358979323846264338327950288419716939937510
typedef double double2[2];
typedef double double3[3];
typedef double double4[4];
typedef double double6[6];
typedef double double9[9];

using namespace cv;
using namespace std;

/*! \brief difference object
*
*   the object which measures the blocks 
*   obtained by the projected graph's difference
*/
struct DIFFObject
{
	int begin;
	int length;
	int area;
	int max_value;
};

/*! \brief defect area
*
*   the object which measures the defect area
*/
typedef struct
{
	int type;
	int size;
	Rect rect;
} DEFECT;


inline int func_nc8(int *b)
/* connectivity detection for each point */
{
	int n_odd[4] = { 1, 3, 5, 7 };  /* odd-number neighbors */
	int i, j, sum, d[10];           /* control variable */

	for (i = 0; i <= 9; i++) 
	{
		j = i;
		if (i == 9) 
			j = 1;
		if (abs(*(b + j)) == 1) 
		{
			d[i] = 1;
		}
		else 
		{
			d[i] = 0;
		}
	}
	sum = 0;
	for (i = 0; i < 4; i++) 
	{
		j = n_odd[i];
		sum = sum + d[j] - d[j] * d[j + 1] * d[j + 2];
	}
	return (sum);
}

inline double get_euclid_dist(Point p, Rect r)
{
	return sqrt(pow(p.x - r.x, 2) + pow(p.y - r.y, 2));
}

inline double get_euclid_dist(double2 p, Rect r)
{
	double2 c;
	c[0] = r.x + r.width / 2;
	c[1] = r.y + r.height / 2;
	return sqrt(pow(p[0] - c[0], 2) + pow(p[1] - c[1], 2));
}

inline double get_euclid_norm(double2 p)
{
	return sqrt(pow(p[0], 2) + pow(p[1], 2));
}


inline double get_euclid_dist(double2 p1, Point p2)
{
	return sqrt(pow(p1[0] - p2.x, 2) + pow(p1[1] - p2.y, 2));
}
inline double get_euclid_dist(Point2f p1, Point2f p2)
{
	return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

inline double get_euclid_dist(Point p1, Point p2)
{
	return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

inline double get_euclid_dist(double* p1, double* p2)
{
	return sqrt(pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2));
}

/*!
segment projection 
according to projection and projection difference
\param projection
\param diff_projection
\param size
\param threshold zero threshold for energy conservation
\param block_set the set of segmented blocks
*/
inline void analysis_diff_projection(int *projection, int *diff_projection
	, int size, int threshold, vector<DIFFObject> &block_set)
{
	int diff_diff;

	int sum = 0;

	DIFFObject object;
	object.area = 0;
	object.begin = 0;
	object.length = 0;
	object.max_value = 0;

	bool begin = true;

	for (int i = 0; i != size+10; ++i)
	{
		if (begin && diff_projection[i] > 0)
		{
			object.begin = i + 1;
			object.area = 0;
			object.length = 0;
			object.max_value = 0;

			begin = false;
		}
		object.area += projection[i];
		object.length++;
		sum+= diff_projection[i];

		if (projection[i] > abs(object.max_value))
		{
			object.max_value = projection[i];
		}

		if (sum < threshold && projection[i]<threshold)
		{
			if (!begin && object.max_value > 10 && object.length > 5)
			{
				block_set.push_back(object);
			}
			
			begin = true;
			sum = 0;
		}
		//else if (i == size - 1)
		//{
		//	block_set.push_back(object);
		//	begin = true;
		//}

	}
}

inline void get_diff_projection(int *projection, int *diff_projection, int size)
{
	for (int i = 0; i != size; ++i)
	{
		diff_projection[i] = projection[(i + 1) % size] - projection[i];
	}
}

/*!
count the number of effective points in ROI.
\param buffer ROI bitmap
\param cols
\param rect
\return the number of effective points
*/
inline int topoStatit(BYTE *buffer,int cols, Rect rect)
{
	int white_num = 0;
	for (int i = rect.x; i != rect.x + rect.width; ++i)
	{
		for (int j = rect.y; j != rect.y + rect.height; ++j)
		{
			if (buffer[j*cols + i] > 200)
			{
				white_num++;
			}
		}
	}

	return white_num;
}

inline int blackPointNum(Mat &mat, Rect rect)
{
	int black_num = 0;
	for (int i = rect.y; i != rect.y + rect.height; ++i)
	{
		uchar *data = (uchar*)mat.ptr<uchar>(i);
		for (int j = rect.x; j != rect.x + rect.width; ++j)
		{
			if (data[j] < 50)
			{
				black_num++;
			}
		}
	}

	return black_num;
}

inline bool isInRect(double2 p, Rect r)
{
	if (p[0] > r.x && p[0]<r.x + r.width && p[1]>r.y && p[1] < r.y + r.height)
	{
		return true;
	}
	else
	{
		return false;
	}
}

inline bool isInRect(Point2f p, Rect r)
{
	if (p.x > r.x && p.x<r.x + r.width && p.y>r.y && p.y < r.y + r.height)
	{
		return true;
	}
	else
	{
		return false;
	}
}

/*!
calculate projection graph, difference graph, 
and then divide projection graph to get character zones.
\param char_area_set character zones
\param buffer gray bitmap
\param cols
\param rect ROI
\param diffThres difference threshold
\param charPointNumThres minimum number of character points
\return number of characters
*/
inline int projectSegment(vector<Rect> &char_area_set,BYTE *buffer,int cols
	, Rect rect,int diffThres,int charPointNumThres)
{
	int projection[10240] = { 0 };
	int diff_projection[10240] = { 0 };

	for (int i = rect.y; i != rect.y + rect.height; ++i)
	{
		for (int j = rect.x; j != rect.x + rect.width; ++j)
		{
			if ( buffer[i*cols+j]> 100)
			{
				projection[j - rect.x + 1]++;
			}
		}
	}

	get_diff_projection(projection, diff_projection, rect.width);

	vector<DIFFObject> block_set;
	analysis_diff_projection(projection, diff_projection, rect.width, diffThres, block_set);

	for (int i = 0; i != block_set.size(); ++i)
	{
		Rect r;
		r.x = block_set[i].begin + rect.x - 2;
		r.y = rect.y - 2;
		r.width = block_set[i].length + 4;
		r.height = rect.height + 4;

		int stat = topoStatit(buffer, cols, r);

		if (charPointNumThres < stat)
		{
			char_area_set.push_back(r);
		}
	}
	return char_area_set.size();
}

//skeleton generation algorithm
inline void HilditchThin(Mat& src, Mat& dst, Rect ROI)
{
	if (src.type() != CV_8UC1)
	{
		printf("Skeletons generating：only for binary image.");
		return;
	}

	if (dst.data != src.data)
	{
		src.copyTo(dst);
	}

	//eight neighborhood offset
	int offset[9][2] = { { 0,0 },{ 1,0 },{ 1,-1 },{ 0,-1 },{ -1,-1 },
	{ -1,0 },{ -1,1 },{ 0,1 },{ 1,1 } };

	int px, py;
	int b[9];                                 
	int counter; 
	int i, x, y,sum;

	int width = dst.cols;
	int height = dst.rows;
	uchar* img = dst.data;
	int step = dst.step;
	do
	{
		counter = 0;

		for (y = ROI.y; y < ROI.y + ROI.height; y++)
		{
			for (x = ROI.x; x < ROI.x + ROI.width; x++)
			{
				for (i = 0; i < 9; i++)
				{
					b[i] = 0;
					px = x + offset[i][0];
					py = y + offset[i][1];
					if (px >= ROI.x 
						&& px < ROI.x + ROI.width 
						&& py >= ROI.y 
						&& py < ROI.y + ROI.height)
					{
						if (img[py*step + px] == WHITE)
						{
							b[i] = 1;
						}
						else if (img[py*step + px] == GRAY)
						{
							b[i] = -1;
						}
					}
				}

				//case 1, is entity
				if (b[0] != 1)
					continue;

				//case 2, is background
				if (abs(b[1]) + abs(b[3]) + abs(b[5]) + abs(b[7]) == 4
					|| abs(b[2]) + abs(b[4]) + abs(b[6]) + abs(b[8]) == 4)
					continue;

				//case 3, endpoint
				if (abs(b[1]) + abs(b[2]) + abs(b[3])
					+ abs(b[4]) + abs(b[5]) + abs(b[6])
					+ abs(b[7]) + abs(b[8])<2)
				{
					continue;
				}

				//case 4, isolated point 
				sum = 0;
				for (i = 1; i <= 8; i++)
				{
					if (b[i] == 1)
					{
						sum++;
					}
				}
				if (sum == 0)
					continue;

				//case 5, continuity
				if (func_nc8(b) != 1)
					continue;

				//case 6, A skeleton of width 2 can only remove 1 edge
				sum = 0;
				for (i = 1; i <= 8; i++)
				{
					if (b[i] != -1)
					{
						sum++;
					}
					else
					{
						b[i] = 0;
						if (func_nc8(b) == 1)
						{
							sum++;
						}

						b[i] = -1;
					}
				}
				if (sum != 8)
					continue;

				img[y*step + x] = GRAY; 
				counter++;
			}
		}

		if (counter != 0)
		{
			for (y = ROI.y; y < ROI.y + ROI.height; y++)
			{
				for (x = ROI.x; x < ROI.x + ROI.width; x++)
				{
					if (img[y*step + x] == GRAY)
					{
						img[y*step + x] = BLACK;
					}
				}
			}
		}

	} while (counter != 0);
}

/*! fill neighborhood of one point */
inline void fillblack(Mat &m, Point2f p, int diaThres)
{
	Point cbTopLeft, cbBottomRight;
	cbTopLeft.x = p.x - diaThres / 2 > 0 ? p.x - diaThres / 2 : 1;
	cbTopLeft.y = p.y - diaThres / 2 > 0 ? p.y - diaThres / 2 : 1;
	cbBottomRight.x = p.x + diaThres / 2 < m.cols
		? p.x + diaThres / 2 : m.cols - 1;

	cbBottomRight.y = p.y + diaThres / 2 < m.rows
		? p.y + diaThres / 2 : m.rows - 1;

	for (int i = cbTopLeft.y; i != cbBottomRight.y; ++i)
	{
		for (int j = cbTopLeft.x; j != cbBottomRight.x; ++j)
		{
			Point2f pinbitmap;
			pinbitmap.x = j;
			pinbitmap.y = i;
			if (diaThres / 2 >= get_euclid_dist(p, pinbitmap))
			{
				if (i < m.rows && j < m.cols)
				{
					m.at<uchar>(i, j) = 0;
				}
			}
		}
	}
}

/*!
calculate projection graph, difference graph,
and then divide projection graph to get character zones.
\param char_area_set character zones
\param buffer gray bitmap
\param cols
\param rect ROI
\param diffThres difference threshold
\param charPointNumThres minimum number of character points
\return number of characters
*/
inline bool colordiff4onepoint(vector<int> &dia_set, int value, double zeroDiaRate, int bI, int eI)
{
	int size = dia_set.size();
	int range = eI - bI;
	int zero_num = 0;
	for (int i = bI; i < eI; ++i)
	{
		if (value < dia_set[i%size])
		{
			continue;
		}

		zero_num++;
	}

	for (int i = bI + size / 2; i < eI + size / 2; ++i)
	{
		if (value < dia_set[i%size])
		{
			continue;
		}

		zero_num++;
	}

	if (zero_num >= range * 2 * zeroDiaRate)
	{
		return true;
	}
	else
	{
		return false;
	}
}

/*!
judging multi-carving by neighbourhood diameter.
\param dia_set a sorted set of diameters
\param diaThres diameter threshold
\return 
*/
inline bool multicarv4onepoint(vector<int> &dia_set, int diaThres)
{
	if (dia_set[0] > diaThres)
	{
		return true;
	}
	else
	{
		return false;
	}
}

inline bool isContain(Rect r_outer, Rect r_inner)
{
	if (r_outer.x <= r_inner.x
		&& r_outer.y <= r_inner.y
		&& r_outer.x + r_outer.width >= r_inner.x + r_inner.width
		&& r_outer.y + r_outer.height >= r_inner.y + r_inner.height)
	{
		return true;
	}
	else
	{
		return false;
	}
}


inline bool SortByM1(const Rect &r1, const Rect &r2)
{
	return r1.x < r2.x;
}

inline void matrix2double3(double3 *d_set, double2 *p_set, int num)
{
	for (int i = 0; i != num; ++i)
	{
		d_set[i][0] = p_set[i][0];
		d_set[i][1] = p_set[i][1];
		d_set[i][2] = 1;
	}
}

inline void double2matrix(double2 *p_set, double3 *d_set, int num)
{
	for (int i = 0; i != num; ++i)
	{
		p_set[i][0] = d_set[i][0];
		p_set[i][1] = d_set[i][1];
	}
}

inline void getboundingbox(Rect &r, vector<Point> &defect_point_set)
{
	r.x = INFINITYNUM;
	r.y = INFINITYNUM;
	int endX = 0;
	int endY = 0;
	for (int i = 0; i != defect_point_set.size(); ++i)
	{
		if (defect_point_set[i].x < r.x)
		{
			r.x = defect_point_set[i].x;
		}

		if (defect_point_set[i].y < r.y)
		{
			r.y = defect_point_set[i].y;
		}

		if (defect_point_set[i].x > endX)
		{
			endX = defect_point_set[i].x;
		}

		if (defect_point_set[i].y > endY)
		{
			endY = defect_point_set[i].y;
		}
	}
	r.width = endX - r.x;
	r.height = endY - r.y;
}

inline void matrix2double3(double3 &d_set, double2 &p_set)
{
	d_set[0] = p_set[0];
	d_set[1] = p_set[1];
	d_set[2] = 1;
}

inline double pointinregion(double2 p, Rect r)
{
	double2 r_center;
	r_center[0] = r.x + r.width / 2;
	r_center[1] = r.y + r.height / 2;
	return get_euclid_dist(p, r_center);
}

template<class TT1>
void  __declspec(dllexport) VectorMinus(TT1 *M, TT1 *M1, TT1 *M2)
{
	M2[0] = M[0] - M1[0];
	M2[1] = M[1] - M1[1];
	M2[2] = M[2] - M1[2];
}


//matrix  M2=M-M1   (double , float) - (double , float) = (double , float)
template<class TT1, class TT2, class TT3>
void  __declspec(dllexport) MatrixMinus(TT1* M, TT2* M1, TT3* M2, int nrows, int ncols)
{//M - nrows*ncols    M1 - nrows*ncols   M2 - nrows*ncols
	for (int i = 0; i < nrows; ++i)
	{
		for (int j = 0; j < ncols; ++j)
		{
			M2[i*ncols + j] = M[i*ncols + j] - M1[i*ncols + j];
		}
	}
}

// (float,double) * (float,double) =  (double , float)
template<class TT1, class TT2, class TT3>
void  __declspec(dllexport) MultMatrix(TT1* M, TT2* M1, TT3* M2, int rows, int cols)
{//M - rows*cols    M1 - cols*rows   M2 - rows*rows
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j<rows; ++j)
		{
			M2[i*rows + j] = 0.0;
			for (int k = 0; k < cols; ++k)
			{
				M2[i*rows + j] += M[i*cols + k] * M1[k*rows + j];
			}
		}
	}
}


// (float,double) * (float,double) =  (double , float)
template<class TT1, class TT2, class TT3>
void  __declspec(dllexport) MultMatrix(TT1* M, TT2* M1, TT3* M2, int rows, int soms, int cols)
{
	//M - rows*soms    M1 - soms*cols   M2 - rows*cols
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j<cols; ++j)
		{
			M2[i*cols + j] = 0.0;
			for (int k = 0; k < soms; ++k)
			{
				M2[i*cols + j] += M[i*soms + k] * M1[k*cols + j];
			}
		}
	}
}


//matrix[3][4] change to  matrix0[4][3]
template<class T, class T0>
void __declspec(dllexport) MatrixTranspose(T* matrix, T0* matrix0, int nrows, int ncols)
{// matrix - nrows*ncols    matrix0 - ncols*nrows
	for (int i = 0; i < nrows; ++i)
	{
		for (int j = 0; j < ncols; ++j)
		{
			matrix0[j*nrows + i] = matrix[i*ncols + j];
		}
	}
}

/// <summary>   
/// Generating identity matrix
/// </summary>     
/// <param name="matrix"></param> 
/// <param name="nrows"></param>
template<class T>
void __declspec(dllexport) MatrixIdentity(T* matrix, int nrows)
{ //  matrix -  nrows*nrows
	memset(matrix, 0, nrows*nrows * sizeof(T));
	for (int i = 0; i < nrows; ++i)
	{
		matrix[i*nrows + i] = 1;
	}
}


/// <summary>   
/// matrix inversion
/// </summary>     
/// <param name="Matrix"></param> 
/// <param name="MatrixA">inverse matrix</param>
/// <param name="rows"></param>
//////////////////////////////////////////////////////////////////////////
// invertible matrice
template<class TT1, class TT2>
void __declspec(dllexport) MatrixAnti(TT1* Matrix, TT2* MatrixA, int rows)
{//  Matrix * MatrixA = I          I = E
	double* E = new double[rows*rows];
	MatrixIdentity(E, rows);

	//Doolittle solution
	Doolittle(Matrix, E, MatrixA, rows);
	delete[]E;
	E = NULL;
}



template<class TT1>
double __declspec(dllexport) getCos(TT1 *v1, TT1 *v2)
{
	double innerProduct = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

	double norm1 = sqrt(pow(v1[0], 2) + pow(v1[1], 2) + pow(v1[2], 2));

	double norm2 = sqrt(pow(v2[0], 2) + pow(v2[1], 2) + pow(v2[2], 2));

	return innerProduct / (norm1*norm2);
}

inline float getCos(Point2f center, Point2f side1, Point2f side2)
{
	Point2f v1, v2;

	v1.x = side1.x - center.x;

	v1.y = side1.y - center.y;

	v2.x = side2.x - center.x;

	v2.y = side2.y - center.y;

	float innerProduct = v1.x * v2.x + v1.y * v2.y;

	float norm1 = sqrt(pow(v1.x, 2) + pow(v1.y, 2));

	float norm2 = sqrt(pow(v2.x, 2) + pow(v2.y, 2));

	return innerProduct / (norm1*norm2);
}

inline float getCos(Point center, Point side1, Point side2)
{
	Point2f v1, v2;

	v1.x = side1.x - center.x;

	v1.y = side1.y - center.y;

	v2.x = side2.x - center.x;

	v2.y = side2.y - center.y;

	float innerProduct = v1.x * v2.x + v1.y * v2.y;

	float norm1 = sqrt(pow(v1.x, 2) + pow(v1.y, 2));

	float norm2 = sqrt(pow(v2.x, 2) + pow(v2.y, 2));

	return innerProduct / (norm1*norm2);
}

template<class TT1>
TT1 __declspec(dllexport) GetSum(TT1 *ps, int ps_num)
{
	TT1 sum = 0;

	for (int i = 0; i != ps_num; ++i)
	{
		sum += ps[i];
	}

	return sum;
}

template<class TT1>
void __declspec(dllexport) linearTransf(TT1(&v)[3], TT1 rotation[9], TT1 transl[3])
{
	TT1 temp[3] = { 0 };
	memcpy(temp, v, sizeof(TT1) * 3);
	MultMatrix(rotation, temp, v, 3, 3, 1);
	v[0] += transl[0];
	v[1] += transl[1];
	v[2] += transl[2];
}

template<class TT1>
void __declspec(dllexport) linearTransf(TT1(*v_set)[3], int N, TT1 rotation[9], TT1 transl[3])
{
	for (int i = 0; i != N; ++i)
	{
		linearTransf(v_set[i], rotation, transl);
	}
}

template<class TT1>
void __declspec(dllexport) linearTransf(TT1(&v)[3], TT1 rotation[9])
{
	TT1 temp[3] = { 0 };
	memcpy(temp, v, sizeof(TT1) * 3);
	MultMatrix(rotation, temp, v, 3, 3, 1);
}

template<class TT1>
void __declspec(dllexport) linearTransf(TT1(*v_set)[3], int N, TT1 rotation[9])
{
	for (int i = 0; i != N; ++i)
	{
		linearTransf(v_set[i], rotation);
	}
}

inline void rotatePoint(Point &p, Point &p_rotate, double degree, int width, int height)
{
	Point2f center;
	center.x = float(width / 2.0 + 0.5);
	center.y = float(height / 2.0 + 0.5);


	double2 point, p1;
	point[0] = p.x - center.x;
	point[1] = p.y - center.y;

	double4 matrix;

	double angle = degree*PI / 180;
	matrix[0] = cos(angle);
	matrix[1] = sin(angle);
	matrix[2] = -1 * sin(angle);
	matrix[3] = cos(angle);

	MultMatrix(point, matrix, p1, 1, 2, 2);

	p_rotate.x = p1[0] + center.x;
	p_rotate.y = p1[1] + center.y;
}

/*!
score points with reference to their relative positions.
\param temp_ctr template center point
\param mea_ctr measured center point
\param thres distance threshold of y component
\return score
*/
inline double ScorePointRela(double2 temp_ctr, double2 mea_ctr, int thres)
{
	double score = 0;

	double distY = abs(temp_ctr[1] - mea_ctr[1]);
	double distX = abs(mea_ctr[0] - temp_ctr[0]);
	if (distY > thres)
	{
		return 0;
	}

	score = 1 / sqrt(pow(distX, 2) + pow(distY, 2));

	return score;
}

/*!
get projection for judging color diff
\param proj projection
\param bar
\param lightThres lower threshold
\param darkThres upper threshold
\param length effective radius
*/
inline void GetProjY4ColorDiff(int *proj, Mat &bar, int lightThres, int darkThres, int length)
{
	memset(proj, 0, sizeof(int)*bar.rows);

	for (int i = 0; i != bar.rows; ++i)
	{
		uchar *data = (uchar*)bar.ptr<uchar>(i);
		for (int j = 0; j != length; ++j)
		{
			if ((data[j] > lightThres && data[j] < darkThres))
			{
				proj[i] ++;
			}
		}
	}
}

/*!
get projection for judging break
\param proj projection
\param bar
\param length effective radius
\param binaryzationType
*/
inline void GetProjY4Break(int *proj, Mat &bar, int length, int binaryzationType)
{
	memset(proj, 0, sizeof(int)*bar.rows);
	if (1 == binaryzationType)
	{
		threshold(bar, bar, 0, 255, THRESH_OTSU | THRESH_BINARY_INV);
	}
	for (int i = 0; i != bar.rows; ++i)
	{
		uchar *data = (uchar*)bar.ptr<uchar>(i);
		for (int j = 0; j != length; ++j)
		{
			if (data[j] > 100)
			{
				proj[i] ++;
			}
		}
	}
}

/*!
get projection for judging defective gum 
\param proj projection
\param bar
*/
inline void GetProjY4Gum(int *proj, Mat &bar)
{
	memset(proj, 0, sizeof(int)*bar.rows);

	for (int i = 0; i != bar.rows; ++i)
	{
		uchar *data = (uchar*)bar.ptr<uchar>(i);
		for (int j = 0; j != bar.cols; ++j)
		{
			if (data[j] > 100)
			{
				proj[i] ++;
			}
			else
			{
				break;
			}
		}
	}
}

/*!
merge rectangular boxes
\param unionSet
\param rect_set
*/
inline void GetUnionSet(Rect &unionSet, vector<Rect> &rect_set)
{
	int left = INFINITYNUM;
	int right = 0;
	int top = INFINITYNUM;
	int bottom = 0;

	for (int i = 0; i != rect_set.size(); ++i)
	{
		if (rect_set[i].x < left)
		{
			left = rect_set[i].x;
		}

		if (rect_set[i].y < top)
		{
			top = rect_set[i].y;
		}

		if (rect_set[i].x + rect_set[i].width > right)
		{
			right = rect_set[i].x + rect_set[i].width;
		}

		if (rect_set[i].y + rect_set[i].height > bottom)
		{
			bottom = rect_set[i].y + rect_set[i].height;
		}
	}

	unionSet.x = left;
	unionSet.y = top;
	unionSet.width = right - left;
	unionSet.height = bottom - top;
}

/*!
calculate the diameter of vertex neighborhood.
\param dia_set diameters sorted from smallest to largest
\param neighbproj
\param size
\param bI begin index
\param eI end index
*/
inline void getdiaset(vector<int> &dia_set, int *neighbproj, int size, int bI, int eI)
{
	dia_set.clear();
	int dia = 0;
	for (int i = bI; i < eI; ++i)
	{
		dia = neighbproj[i] + neighbproj[i + size / 2];
		dia_set.push_back(dia);
	}

	sort(dia_set.begin(), dia_set.end());
}

/*!
calculate the radius of vertex neighborhood.
\param dia_set radius set
\param neighbproj
*/
inline void getradiusset4break(vector<int> &radius_set, int *neighbproj)
{
	radius_set.clear();
	int dia = 0;
	for (int i = 0; i != 360; ++i)
	{
		dia = neighbproj[i];
		radius_set.push_back(dia);
	}

	//sort(dia_set.begin(), dia_set.end());
}

static bool SortByX4Asc(const Rect &v1, const Rect &v2)
{
	return v1.x < v2.x;//ascending sort
}

static bool SortByHeightDesc(const Rect &r1, const Rect &r2)
{
	return r1.height > r2.height;
}

static bool SortByAreaDesc(const Rect &r1, const Rect &r2)
{
	return r1.height*r1.width > r2.height*r2.width;
}

/*!
filter SN code through boundingbox,
sort by area and orientation.
\param area_set
\param boundingbox
\param r_set
\param num
\param inner_num the number of sn code
\return 
*/
inline bool GetAreaSetFromBoundingBox(vector<Rect> &area_set, Rect &boundingbox
	, Rect *r_set, int num, int inner_num)
{
	for (int i = 0; i != num; ++i)
	{
		if (isContain(boundingbox, r_set[i]))
		{
			area_set.push_back(r_set[i]);
		}
	}
	int index = 0;
	if (inner_num == area_set.size())
	{
		sort(area_set.begin(), area_set.end(), SortByX4Asc);
		return true;
	}
	else if (inner_num < area_set.size())
	{
		index = area_set.size() - inner_num;
		sort(area_set.begin(), area_set.end(), SortByAreaDesc);
		sort(area_set.begin(), area_set.end() - index, SortByX4Asc);
		return true;
	}
	else
	{
		return false;
	}
}

/*!
get boundingbox of point set 
\param boundingbox
\param point_set
\param num
\param margin
*/
inline void getboundingbox(Rect &boundingbox, double2 *point_set, int num, int margin)
{
	boundingbox.x = INFINITYNUM;
	boundingbox.y = INFINITYNUM;

	double2 bottomRight;
	bottomRight[0] = 0;
	bottomRight[1] = 0;

	for (int i = 0; i != num; ++i)
	{
		if (boundingbox.x > point_set[i][0])
		{
			boundingbox.x = point_set[i][0];
		}

		if (boundingbox.y > point_set[i][1])
		{
			boundingbox.y = point_set[i][1];
		}

		if (bottomRight[0] < point_set[i][0])
		{
			bottomRight[0] = point_set[i][0];
		}

		if (bottomRight[1] < point_set[i][1])
		{
			bottomRight[1] = point_set[i][1];
		}
	}
	boundingbox.x -= margin;
	boundingbox.y -= margin;
	boundingbox.width = bottomRight[0] - boundingbox.x + margin;
	boundingbox.height = bottomRight[1] - boundingbox.y + margin;
}

/*!
get boundingbox of point set
\param topLeft
\param bottomRight
\param ps
\param num
\return boundingbox
*/
inline Rect getrect(Point &topLeft, Point &bottomRight, Point *ps, int num)
{
	Rect rect;
	topLeft.x = INFINITYNUM;
	topLeft.y = INFINITYNUM;
	bottomRight.x = 0;
	bottomRight.y = 0;
	for (int i = 0; i != num; ++i)
	{
		if (ps[i].x < topLeft.x)
		{
			topLeft.x = ps[i].x;
		}

		if (ps[i].y < topLeft.y)
		{
			topLeft.y = ps[i].y;
		}

		if (ps[i].x > bottomRight.x)
		{
			bottomRight.x = ps[i].x;
		}

		if (ps[i].y > bottomRight.y)
		{
			bottomRight.y = ps[i].y;
		}
	}

	rect.x = topLeft.x;
	rect.y = topLeft.y;
	rect.width = bottomRight.x - topLeft.x;
	rect.height = bottomRight.y - topLeft.y;

	return rect;
}

/*! classify these ng points into ok and ng according to eliminating outliers */
inline bool elimioutClassify(vector<Point> &ps, int radius, int ngCountThres)
{
	vector<Point>::iterator iter;
	vector<Point>::iterator iter_inner;
	int innerCount = 0;
	for (iter = ps.begin(); iter != ps.end(); ++iter)
	{
		innerCount = 0;
		for (iter_inner = ps.begin(); iter_inner != ps.end(); ++iter_inner)
		{
			if (iter == iter_inner)
			{
				continue;
			}

			if (radius > get_euclid_dist(*iter, *iter_inner))
			{
				innerCount++;
			}

			if (ngCountThres <= innerCount)
			{
				return true;
			}
		}
	}

	return false;
}
#endif
