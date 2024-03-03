/*
Copyright 2011. All rights reserved.
Institute of Measurement and Control Systems
Karlsruhe Institute of Technology, Germany

Authors: Andreas Geiger

libicp is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or any later version.

libicp is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
libicp; if not, write to the Free Software Foundation, Inc., 51 Franklin
Street, Fifth Floor, Boston, MA 02110-1301, USA 
*/

#include "icp.h"

using namespace std;

Icp::Icp(CPoint2 *&ps, int &p_num) :
	dim(2), sub_step(2), max_iter(200), min_delta(1e-6)
{
	double *M = new double[p_num * 2];
	int M_num = p_num;

	pointset2T(M, ps, p_num);
	// check for minimum number of points
	if (M_num < 5) {
		//cout << "ERROR: LIBICP works only with at least 5 model points" << endl;
		M_tree = 0;
		return;
	}

	// copy model points to M_data
	M_data.resize(boost::extents[M_num][dim]);
	for (int m = 0; m < M_num; m++)
		for (int n = 0; n < dim; n++)
			M_data[m][n] = (float)M[m*dim + n];

	// build a kd tree from the model point cloud
	M_tree = new kdtree::KDTree(M_data);

	delete[]M;
}

Icp::Icp(CPoint3 *&ps, int &p_num) :
dim(3), sub_step(2), max_iter(200), min_delta(1e-6)
{
	double *M = new double[p_num * 3];
	int M_num = p_num;

	pointset3T(M, ps, p_num);
	// check for minimum number of points
	if (M_num < 5) {
		//cout << "ERROR: LIBICP works only with at least 5 model points" << endl;
		M_tree = 0;
		return;
	}

	// copy model points to M_data
	M_data.resize(boost::extents[M_num][dim]);
	for (int m = 0; m < M_num; m++)
		for (int n = 0; n < dim; n++)
			M_data[m][n] = (float)M[m*dim + n];

	// build a kd tree from the model point cloud
	M_tree = new kdtree::KDTree(M_data);

	delete[]M;
}

Icp::Icp(CPoint2 *&ps, int &p_num, int iterCount) :
	dim(2), sub_step(2), max_iter(iterCount), min_delta(1e-6)
{
	double *M = new double[p_num * 2];
	int M_num = p_num;

	pointset2T(M, ps, p_num);
	// check for minimum number of points
	if (M_num < 5) {
		//cout << "ERROR: LIBICP works only with at least 5 model points" << endl;
		M_tree = 0;
		return;
	}

	// copy model points to M_data
	M_data.resize(boost::extents[M_num][dim]);
	for (int m = 0; m < M_num; m++)
		for (int n = 0; n < dim; n++)
			M_data[m][n] = (float)M[m*dim + n];

	// build a kd tree from the model point cloud
	M_tree = new kdtree::KDTree(M_data);

	delete[]M;
}


Icp::Icp(CPoint3 *&ps, int &p_num, int iterCount) :
dim(3), sub_step(2), max_iter(iterCount), min_delta(1e-6)
{
	double *M = new double[p_num * 3];
	int M_num = p_num;

	pointset3T(M, ps, p_num);
	// check for minimum number of points
	if (M_num < 5) {
		//cout << "ERROR: LIBICP works only with at least 5 model points" << endl;
		M_tree = 0;
		return;
	}

	// copy model points to M_data
	M_data.resize(boost::extents[M_num][dim]);
	for (int m = 0; m < M_num; m++)
		for (int n = 0; n < dim; n++)
			M_data[m][n] = (float)M[m*dim + n];

	// build a kd tree from the model point cloud
	M_tree = new kdtree::KDTree(M_data);

	delete[]M;
}


Icp::Icp(double *M, const int M_num, const int dim) :
	dim(dim), sub_step(10), max_iter(200), min_delta(1e-6) 
{

	// check for correct dimensionality
	if (dim != 2 && dim != 3) {
		//cout << "ERROR: LIBICP works only for data of dimensionality 2 or 3" << endl;
		M_tree = 0;
		return;
	}

	// check for minimum number of points
	if (M_num < 5) {
		//cout << "ERROR: LIBICP works only with at least 5 model points" << endl;
		M_tree = 0;
		return;
	}

	// copy model points to M_data
	M_data.resize(boost::extents[M_num][dim]);
	for (int m = 0; m < M_num; m++)
		for (int n = 0; n < dim; n++)
			M_data[m][n] = (float)M[m*dim + n];

	// build a kd tree from the model point cloud
	M_tree = new kdtree::KDTree(M_data);
}

Icp::~Icp() {
	if (M_tree)
		delete M_tree;
}

void Icp::pointset2T(double *T, CPoint2 *ps, const int p_num)
{
	for (int i = 0; i != p_num; ++i)
	{
		T[i * 2] = ps[i][0];
		T[i * 2 + 1] = ps[i][1];
	}
}

void Icp::pointset3T(double *T, CPoint3 *ps, const int p_num)
{
	for (int i = 0; i != p_num; ++i)
	{
		T[i * 3] = ps[i][0];
		T[i * 3 + 1] = ps[i][1];
		T[i * 3 + 2] = ps[i][2];
	}
}

void Icp::fit(CPoint2 *&ps, int &p_num, Matrix &R, Matrix &t, vector<int> &fit_set, const double indist)
{
	double *T = new double[p_num * 2];
	int T_num = p_num;

	pointset2T(T, ps, p_num);

	fit(T, T_num, R, t, fit_set, indist);

	delete[] T;
}

void Icp::fit(CPoint3 *&ps, int &p_num, Matrix &R, Matrix &t, vector<int> &fit_set, const double indist)
{
	double *T = new double[p_num * 3];
	int T_num = p_num;

	pointset3T(T, ps, p_num);

	fit(T, T_num, R, t, fit_set, indist);

	delete[] T;
}

void Icp::fit(double *T, const int T_num, Matrix &R, Matrix &t, vector<int> &fit_set, const double indist) {

	// make sure we have a model tree
	if (!M_tree) 
	{
		//cout << "ERROR: No model available." << endl;
		return;
	}

	// check for minimum number of points
	if (T_num < 5) 
	{
		//cout << "ERROR: Icp works only with at least 5 template points" << endl;
		return;
	}
	vector<int> active;
	for (int i = 0; i < T_num; i++)
	{
		active.push_back(i);
	}

	fitIterate(T, T_num, R, t, active, fit_set);
	//getInliers(fit_set, T, T_num, R, t, indist);
}

void Icp::fitIterate(double *&T, const int &T_num, Matrix &R, Matrix &t, const vector<int> active,vector<int> &fit_set)
{
	// check if we have at least 5 active points
	if (active.size() < 5)
		return;

	// iterate until convergence
	for (int iter = 0; iter < max_iter; iter++)
	{
		if (fitStep(T, T_num, R, t, active,fit_set) < min_delta)
		{
			break;
		}
	}
}
