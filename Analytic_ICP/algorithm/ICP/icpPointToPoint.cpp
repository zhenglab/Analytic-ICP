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

#include "icpPointToPoint.h"

using namespace std;

// Also see (3d part): "Least-Squares Fitting of Two 3-D Point Sets" (Arun, Huang and Blostein)
double IcpPointToPoint::fitStep(double *T, const int T_num, MatrixXd &R, VectorXd &t, const vector<int> active, vector<int> &fit_set)
{
	int active_num = active.size();
  // kd tree query + result
  std::vector<float>         query(dim);
  kdtree::KDTreeResultVector result;
  
  // init matrix for point correspondences
  MatrixXd p_m = MatrixXd::Zero(active_num, dim); // model
  MatrixXd p_t = MatrixXd::Zero(active_num, dim); // template
  
  // init mean
  MatrixXd mu_m = MatrixXd::Zero(1, dim);
  MatrixXd mu_t = MatrixXd::Zero(1, dim);
  fit_set.clear();
  // dimensionality 2
  if (dim == 2) {

	  // extract matrix and translation vector
	  double r00 = R(0, 0); double r01 = R(0, 1);
	  double r10 = R(1, 0); double r11 = R(1, 1);
	  double t0 = t(0); double t1 = t(1);

	  // establish correspondences
	  for (int i = 0; i < active_num; i++)
	  {
		  // get index of active point
		  int idx = active[i];
		  
		  // transform point according to R|t
		  query[0] = r00*T[idx * 2 + 0] + r01*T[idx * 2 + 1] + t0;
		  query[1] = r10*T[idx * 2 + 0] + r11*T[idx * 2 + 1] + t1;

		  // search nearest neighbor
		  M_tree->n_nearest(query, 1, result);
		  fit_set.push_back(result[0].idx);
		  // set model point
		  p_m(i, 0) = M_tree->the_data[result[0].idx][0]; mu_m(0,0) += p_m(i, 0);
		  p_m(i, 1) = M_tree->the_data[result[0].idx][1]; mu_m(0,1) += p_m(i, 1);

		  // set template point
		  p_t(i, 0) = query[0]; mu_t(0, 0) += p_t(i, 0);
		  p_t(i, 1) = query[1]; mu_t(0, 1) += p_t(i, 1);
	  }

	  // dimensionality 3
  }
  else
  {

	  // extract matrix and translation vector
	  double r00 = R(0,0); double r01 = R(0,1); double r02 = R(0,2);
	  double r10 = R(1,0); double r11 = R(1,1); double r12 = R(1,2);
	  double r20 = R(2,0); double r21 = R(2,1); double r22 = R(2,2);
	  double t0 = t(0); double t1 = t(1); double t2 = t(2);

	  // establish correspondences
	  for (int i = 0; i < active_num; i++) {

		  // get index of active point
		  int idx = active[i];
		  // transform point according to R|t
		  query[0] = r00*T[idx * 3 + 0] + r01*T[idx * 3 + 1] + r02*T[idx * 3 + 2] + t0;
		  query[1] = r10*T[idx * 3 + 0] + r11*T[idx * 3 + 1] + r12*T[idx * 3 + 2] + t1;
		  query[2] = r20*T[idx * 3 + 0] + r21*T[idx * 3 + 1] + r22*T[idx * 3 + 2] + t2;

		  // search nearest neighbor
		  M_tree->n_nearest(query, 1, result);
		  fit_set.push_back(result[0].idx);
		  // set model point
		  p_m(i, 0) = M_tree->the_data[result[0].idx][0]; mu_m(0) += p_m(i, 0);
		  p_m(i, 1) = M_tree->the_data[result[0].idx][1]; mu_m(1) += p_m(i, 1);
		  p_m(i, 2) = M_tree->the_data[result[0].idx][2]; mu_m(2) += p_m(i, 2);

		  // set template point
		  p_t(i, 0) = query[0]; mu_t(0) += p_t(i, 0);
		  p_t(i, 1) = query[1]; mu_t(1) += p_t(i, 1);
		  p_t(i, 2) = query[2]; mu_t(2) += p_t(i, 2);
	  }
  }
  
  // subtract mean
  mu_m = mu_m/(double)active_num;
  mu_t = mu_t/(double)active_num;
  MatrixXd q_m = p_m - MatrixXd::Ones(active_num, 1)*mu_m;
  MatrixXd q_t = p_t - MatrixXd::Ones(active_num, 1)*mu_t;

  // compute relative rotation matrix R and translation vector t
  MatrixXd H = q_t.transpose()*q_m;
  MatrixXd U, W, V;
  JacobiSVD<MatrixXd> svd(H, ComputeThinU | ComputeThinV);

  U = svd.matrixU();

  V = svd.matrixV();

  W = svd.singularValues();
  MatrixXd R_ = V*U.transpose();
  MatrixXd t_ = mu_m.transpose() - R_*mu_t.transpose();
  
  // compose: R|t = R_|t_ * R|t
  R = R_*R;
  t = R_*t+t_;

  // return max delta in parameters
  if (dim==2) 
	  return max((R_ - MatrixXd::Identity(2,2)).norm(), t_.norm());
  else        
	  return max((R_ - MatrixXd::Identity(3, 3)).norm(), t_.norm());
}

void  IcpPointToPoint::getInliers(vector<int> &fit_set, double *T, const int T_num, const MatrixXd &R, const VectorXd &t, const double indist)
{
	// init inlier vector + query point + query result
	vector<int>            inliers;
	std::vector<float>         query(dim);
	kdtree::KDTreeResultVector neighbor;
	fit_set.clear();

	// dimensionality 2
	if (dim == 2) 
	{

		// extract matrix and translation vector
		double r00 = R(0, 0); double r01 = R(0, 1);
		double r10 = R(1, 0); double r11 = R(1, 1);
		double t0 = t(0); double t1 = t(1);

		// check for all points if they are inliers
		for (int i = 0; i < T_num; i++)
		{

			// transform point according to R|t
			query[0] = r00*T[i * 2 + 0] + r01*T[i * 2 + 1] + t0;
			query[1] = r10*T[i * 2 + 0] + r11*T[i * 2 + 1] + t1;

			// search nearest neighbor
			M_tree->n_nearest(query, 1, neighbor);

			// check if it is an inlier
			if (neighbor[0].dis < indist)
			{
				fit_set.push_back(neighbor[0].idx);
			}
		}

		// dimensionality 3
	}
	else 
	{

		// extract matrix and translation vector
		double r00 = R(0, 0); double r01 = R(0, 1); double r02 = R(0, 2);
		double r10 = R(1, 0); double r11 = R(1, 1); double r12 = R(1, 2);
		double r20 = R(2, 0); double r21 = R(2, 1); double r22 = R(2, 2);
		double t0 = t(0); double t1 = t(1); double t2 = t(2);

		// check for all points if they are inliers
		for (int i = 0; i < T_num; i++)
		{

			// transform point according to R|t
			query[0] = r00*T[i * 3 + 0] + r01*T[i * 3 + 1] + r02*T[i * 3 + 2] + t0;
			query[1] = r10*T[i * 3 + 0] + r11*T[i * 3 + 1] + r12*T[i * 3 + 2] + t1;
			query[2] = r20*T[i * 3 + 0] + r21*T[i * 3 + 1] + r22*T[i * 3 + 2] + t2;

			// search nearest neighbor
			M_tree->n_nearest(query, 1, neighbor);

			// check if it is an inlier
			if (neighbor[0].dis < indist)
			{
				fit_set.push_back(neighbor[0].idx);
			}
		}
	}
}
