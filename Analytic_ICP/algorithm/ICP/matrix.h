/*
Copyright 2011. All rights reserved.
Institute of Measurement and Control Systems
Karlsruhe Institute of Technology, Germany

Authors: Andreas Geiger

matrix is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or any later version.

matrix is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
matrix; if not, write to the Free Software Foundation, Inc., 51 Franklin
Street, Fifth Floor, Boston, MA 02110-1301, USA 
*/

#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <assert.h>

#define endll endl << endl // double end line definition

class Matrix {

public:

  // constructor / deconstructor
  Matrix ();                                                  // init empty 0x0 matrix
  Matrix (const int m,const int n);                   // init empty mxn matrix
  Matrix (const int m,const int n,const double* val_); // init mxn matrix with values from array 'val'
  Matrix (const Matrix &M);                                   // creates deepcopy of M
  ~Matrix ();

  // assignment operator, copies contents of M
  Matrix& operator= (const Matrix &M);

  // copies submatrix of M into array 'val', default values copy whole row/column/matrix
  void getData(double* val_, int i1=0, int j1=0, int i2=-1, int j2=-1);

  // set or get submatrices of current matrix
  Matrix getMat(int i1, int j1, int i2=-1, int j2=-1);
  void   setMat(const Matrix &M,const int i,const int j);

  // set sub-matrix to scalar (default 0), -1 as end replaces whole row/column/matrix
  void setVal(double s, int i1=0, int j1=0, int i2=-1, int j2=-1);

  // set (part of) diagonal to scalar, -1 as end replaces whole diagonal
  void setDiag(double s, int i1=0, int i2=-1);

  // clear matrix
  void zero();
  
  // extract columns with given index
  Matrix extractCols (std::vector<int> idx);

  // create identity matrix
  static Matrix eye (const int m);
  void          eye ();

  // create matrix with ones
  static Matrix ones(const int m,const int n);

  // create diagonal matrix with nx1 or 1xn matrix M as elements
  static Matrix diag(const Matrix &M);
  
  // returns the m-by-n matrix whose elements are taken column-wise from M
  static Matrix reshape(const Matrix &M, int m, int n);

  // create 3x3 rotation matrices (convention: http://en.wikipedia.org/wiki/Rotation_matrix)
  static Matrix rotMatX(const double &angle);
  static Matrix rotMatY(const double &angle);
  static Matrix rotMatZ(const double &angle);

  // simple arithmetic operations
  Matrix  operator+ (const Matrix &M); // add matrix
  Matrix  operator- (const Matrix &M); // subtract matrix
  Matrix  operator* (const Matrix &M); // multiply with matrix
  Matrix  operator* (const double &s);  // multiply with scalar
  Matrix  operator/ (const Matrix &M); // divide elementwise by matrix (or vector)
  Matrix  operator/ (const double &s);  // divide by scalar
  Matrix  operator- ();                // negative matrix
  Matrix  operator~ ();                // transpose
  double   l2norm ();                   // euclidean norm (vectors) / frobenius norm (matrices)
  double   mean ();                     // mean of all elements in matrix

  // complex arithmetic operations
  static Matrix cross (const Matrix &a, const Matrix &b);    // cross product of two vectors
  static Matrix inv (const Matrix &M);                       // invert matrix M
  bool   inv ();                                             // invert this matrix
  double  det ();                                             // returns determinant of matrix
  bool   solve (const Matrix &M,double eps=1e-20);            // solve linear system M*x=B, replaces *this and M
  bool   lu(int *idx, double &d, double eps=1e-20);        // replace *this by lower upper decomposition
  void   svd(Matrix &U,Matrix &W,Matrix &V);                 // singular value decomposition *this = U*diag(W)*V^T

  // print matrix to stream
  friend std::ostream& operator<< (std::ostream& out,const Matrix& M);

  // direct data access
  double   **val;
  int   m,n;

private:

  void allocateMemory (const int m_,const int n_);
  void releaseMemory ();
  inline double pythag(double a,double b);

};

class CPoint2 :public Matrix
{
public:
	CPoint2() :Matrix(2, 1)
	{
		val[0][0] = 0;
		val[1][0] = 0;
	}
	/*!
	*	reference to the values, CPoint2[0] x value ; CPoint2[1] y value;
	*/
	double & operator[](int i) { assert(i < 2 && i >= 0); return val[i][0]; };

	/*!
	*	constant values, CPoint2[0] x value ; CPoint2[1] y value;
	*/
	double  operator[](int i) const { assert(i < 2 && i >= 0); return val[i][0]; };

	CPoint2 operator* (const Matrix &M) 
	{
		const CPoint2 &A = *this;
		const Matrix &B = M;
		if (A.n != B.m) 
		{
			printf("ERROR: Trying to multiply matrices of size (%fx%f) and (%fx%f)", A.m, A.n, B.m, B.n);
			exit(0);
		}
		CPoint2 C;
		for (int i = 0; i<A.m; i++)
			for (int j = 0; j<B.n; j++)
				for (int k = 0; k<A.n; k++)
					C.val[i][j] += A.val[i][k] * B.val[k][j];
		return C;
	}

	CPoint2 operator+ (const Matrix &M) {
		const CPoint2 &A = *this;
		const Matrix &B = M;
		if (A.m != B.m || A.n != B.n) 
		{
			printf("ERROR: Trying to add matrices of size (%fx%f) and (%fx%f)", A.m, A.n, B.m, B.n);
			exit(0);
		}
		CPoint2 C;
		for (int i = 0; i<m; i++)
			for (int j = 0; j<n; j++)
				C.val[i][j] = A.val[i][j] + B.val[i][j];
		return C;
	}

	CPoint2 operator- (const CPoint2 &M) {
		const CPoint2 &A = *this;
		const CPoint2 &B = M;
		CPoint2 C;
		C.val[0][0] = A.val[0][0] - B.val[0][0];
		C.val[1][0] = A.val[1][0] - B.val[1][0];
		return C;
	}

	CPoint2& operator=(const CPoint2& cls)
	{
		val[0][0] = cls.val[0][0];
		val[1][0] = cls.val[1][0];

		return *this;
	}

	CPoint2& operator=(const Matrix& cls)
	{
		if (cls.m != 2 || cls.n != 1)
		{
			printf("ERROR: Trying to assign matrices of size (%fx%f)", cls.m, cls.n);
			exit(0);
		}
		val[0][0] = cls.val[0][0];
		val[1][0] = cls.val[1][0];

		return *this;
	}

};


class CPoint3 :public Matrix
{
public:
	CPoint3() :Matrix(3, 1)
	{
		val[0][0] = 0;
		val[1][0] = 0;
		val[2][0] = 0;
	}
	/*!
	*	reference to the values, CPoint2[0] x value ; CPoint2[1] y value;
	*/
	double & operator[](int i) { assert(i < 3 && i >= 0); return val[i][0]; };

	/*!
	*	constant values, CPoint2[0] x value ; CPoint2[1] y value;
	*/
	double  operator[](int i) const { assert(i < 3 && i >= 0); return val[i][0]; };

	CPoint3 operator* (const Matrix &M)
	{
		const CPoint3 &A = *this;
		const Matrix &B = M;
		if (A.n != B.m)
		{
			printf("ERROR: Trying to multiply matrices of size (%fx%f) and (%fx%f)"
				, A.m, A.n, B.m, B.n);
			exit(0);
		}
		CPoint3 C;
		for (int i = 0; i<A.m; i++)
			for (int j = 0; j<B.n; j++)
				for (int k = 0; k<A.n; k++)
					C.val[i][j] += A.val[i][k] * B.val[k][j];
		return C;
	}

	CPoint3 operator+ (const Matrix &M) {
		const CPoint3 &A = *this;
		const Matrix &B = M;
		if (A.m != B.m || A.n != B.n)
		{
			printf("ERROR: Trying to add matrices of size (%fx%f) and (%fx%f)"
				, A.m, A.n, B.m, B.n);
			exit(0);
		}
		CPoint3 C;
		for (int i = 0; i<m; i++)
			for (int j = 0; j<n; j++)
				C.val[i][j] = A.val[i][j] + B.val[i][j];
		return C;
	}

	CPoint3 operator- (const CPoint3 &M) {
		const CPoint3 &A = *this;
		const CPoint3 &B = M;
		CPoint3 C;
		C.val[0][0] = A.val[0][0] - B.val[0][0];
		C.val[1][0] = A.val[1][0] - B.val[1][0];
		C.val[2][0] = A.val[2][0] - B.val[2][0];
		return C;
	}

	CPoint3& operator=(const CPoint3& cls)
	{
		val[0][0] = cls.val[0][0];
		val[1][0] = cls.val[1][0];
		val[2][0] = cls.val[2][0];

		return *this;
	}

	CPoint3& operator=(const Matrix& cls)
	{
		if (cls.m != 3 || cls.n != 1)
		{
			printf("ERROR: Trying to assign matrices of size (%fx%f)", cls.m, cls.n);
			exit(0);
		}
		val[0][0] = cls.val[0][0];
		val[1][0] = cls.val[1][0];
		val[2][0] = cls.val[2][0];

		return *this;
	}

};
#endif // MATRIX_H
