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

#include "icpPointToPlane.h"

using namespace std;

// Also see (3d part): "Linear Least-Squares Optimization for Point-to-Plane ICP Surface Registration" (Kok-Lim Low)
double IcpPointToPlane::fitStep(double *T, const int T_num, Matrix &R, Matrix &t, const vector<int> active, vector<int> &fit_set)
{
	int active_num = active.size();
	// kd tree query + result
	std::vector<float>         query(dim);
	kdtree::KDTreeResultVector result;

	// init matrix for point correspondences
	Matrix p_m(active_num, dim); // model
	Matrix p_t(active_num, dim); // template
	fit_set.clear();
	// dimensionality 2
	if (dim == 2) {

		// extract matrix and translation vector
		double r00 = R.val[0][0]; double r01 = R.val[0][1];
		double r10 = R.val[1][0]; double r11 = R.val[1][1];
		double t0 = t.val[0][0]; double t1 = t.val[1][0];

		// init A and b
		Matrix A(active_num, 3);
		Matrix b(active_num, 1);

		// establish correspondences
		for (int i = 0; i < active_num; i++) {

			// get index of active point
			int idx = active[i];
			// transform point according to R|t
			query[0] = r00*T[idx * 2 + 0] + r01*T[idx * 2 + 1] + t0;
			query[1] = r10*T[idx * 2 + 0] + r11*T[idx * 2 + 1] + t1;

			// search nearest neighbor
			M_tree->n_nearest(query, 1, result);
			fit_set.push_back(result[0].idx);
			// model point
			double dx = M_tree->the_data[result[0].idx][0];
			double dy = M_tree->the_data[result[0].idx][1];

			// model point normal
			double nx = M_normal[result[0].idx * 2 + 0];
			double ny = M_normal[result[0].idx * 2 + 1];

			// template point
			double sx = query[0];
			double sy = query[1];

			// setup least squares system
			A.val[i][0] = ny*sx - nx*sy;
			A.val[i][1] = nx;
			A.val[i][2] = ny;
			b.val[i][0] = nx*dx + ny*dy - nx*sx - ny*sy;
		}

		// linear least square matrices
		Matrix A_ = ~A*A;
		Matrix b_ = ~A*b;

		// solve linear system
		if (b_.solve(A_)) {

			// rotation matrix
			Matrix R_ = Matrix::eye(2);
			R_.val[0][1] = -b_.val[0][0];
			R_.val[1][0] = +b_.val[0][0];

			// orthonormalized rotation matrix
			Matrix U, W, V;
			R_.svd(U, W, V);
			R_ = U*~V;

			// translation vector
			Matrix t_(2, 1);
			t_.val[0][0] = b_.val[1][0];
			t_.val[1][0] = b_.val[2][0];

			// compose: R|t = R_|t_ * R|t
			R = R_*R;
			t = R_*t + t_;
			return max((R_ - Matrix::eye(2)).l2norm(), t_.l2norm());
		}

		// dimensionality 3
	}
	else {

		// extract matrix and translation vector
		double r00 = R.val[0][0]; double r01 = R.val[0][1]; double r02 = R.val[0][2];
		double r10 = R.val[1][0]; double r11 = R.val[1][1]; double r12 = R.val[1][2];
		double r20 = R.val[2][0]; double r21 = R.val[2][1]; double r22 = R.val[2][2];
		double t0 = t.val[0][0]; double t1 = t.val[1][0]; double t2 = t.val[2][0];

		// init A and b
		Matrix A(active_num, 6);
		Matrix b(active_num, 1);

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
			// model point
			double dx = M_tree->the_data[result[0].idx][0];
			double dy = M_tree->the_data[result[0].idx][1];
			double dz = M_tree->the_data[result[0].idx][2];

			// model point normal
			double nx = M_normal[result[0].idx * 3 + 0];
			double ny = M_normal[result[0].idx * 3 + 1];
			double nz = M_normal[result[0].idx * 3 + 2];

			// template point
			double sx = query[0];
			double sy = query[1];
			double sz = query[2];

			// setup least squares system
			A.val[i][0] = nz*sy - ny*sz;
			A.val[i][1] = nx*sz - nz*sx;
			A.val[i][2] = ny*sx - nx*sy;
			A.val[i][3] = nx;
			A.val[i][4] = ny;
			A.val[i][5] = nz;
			b.val[i][0] = nx*dx + ny*dy + nz*dz - nx*sx - ny*sy - nz*sz;
		}

		// linear least square matrices
		Matrix A_ = ~A*A;
		Matrix b_ = ~A*b;

		// solve linear system
		if (b_.solve(A_)) {

			// rotation matrix
			Matrix R_ = Matrix::eye(3);
			R_.val[0][1] = -b_.val[2][0];
			R_.val[1][0] = +b_.val[2][0];
			R_.val[0][2] = +b_.val[1][0];
			R_.val[2][0] = -b_.val[1][0];
			R_.val[1][2] = -b_.val[0][0];
			R_.val[2][1] = +b_.val[0][0];

			// orthonormalized rotation matrix
			Matrix U, W, V;
			R_.svd(U, W, V);
			R_ = U*~V;

			// translation vector
			Matrix t_(3, 1);
			t_.val[0][0] = b_.val[3][0];
			t_.val[1][0] = b_.val[4][0];
			t_.val[2][0] = b_.val[5][0];

			// compose: R|t = R_|t_ * R|t
			R = R_*R;
			t = R_*t + t_;
			return max((R_ - Matrix::eye(3)).l2norm(), t_.l2norm());
		}
	}

	// failure
	return 0;
}

void IcpPointToPlane::getInliers(vector<int> &fit_set, double *T, const int T_num, const Matrix &R, const Matrix &t, const double indist) {

	// init inlier vector + query point + query result
	vector<int>            inliers;
	std::vector<float>         query(dim);
	kdtree::KDTreeResultVector neighbor;
	fit_set.clear();

	// dimensionality 2
	if (dim == 2) {

		// extract matrix and translation vector
		double r00 = R.val[0][0]; double r01 = R.val[0][1];
		double r10 = R.val[1][0]; double r11 = R.val[1][1];
		double t0 = t.val[0][0]; double t1 = t.val[1][0];

		// check for all points if they are inliers
		for (int i = 0; i < T_num; i++) {

			// transform point according to R|t
			double sx = r00*T[i * 2 + 0] + r01*T[i * 2 + 1]; query[0] = sx;
			double sy = r10*T[i * 2 + 0] + r11*T[i * 2 + 1]; query[1] = sy;

			// search nearest neighbor
			M_tree->n_nearest(query, 1, neighbor);

			// model point
			double dx = M_tree->the_data[neighbor[0].idx][0];
			double dy = M_tree->the_data[neighbor[0].idx][1];

			// model point normal
			double nx = M_normal[neighbor[0].idx * 2 + 0];
			double ny = M_normal[neighbor[0].idx * 2 + 1];

			// check if it is an inlier
			if ((sx - dx)*nx + (sy - dy)*ny < indist)
			{
				fit_set.push_back(neighbor[0].idx);
			}
		}

		// dimensionality 3
	}
	else {

		// extract matrix and translation vector
		double r00 = R.val[0][0]; double r01 = R.val[0][1]; double r02 = R.val[0][2];
		double r10 = R.val[1][0]; double r11 = R.val[1][1]; double r12 = R.val[1][2];
		double r20 = R.val[2][0]; double r21 = R.val[2][1]; double r22 = R.val[2][2];
		double t0 = t.val[0][0]; double t1 = t.val[1][0]; double t2 = t.val[2][0];

		// check for all points if they are inliers
		for (int i = 0; i < T_num; i++) {

			// transform point according to R|t
			double sx = r00*T[i * 3 + 0] + r01*T[i * 3 + 1] + r02*T[i * 3 + 2] + t0; query[0] = sx;
			double sy = r10*T[i * 3 + 0] + r11*T[i * 3 + 1] + r12*T[i * 3 + 2] + t1; query[1] = sy;
			double sz = r20*T[i * 3 + 0] + r21*T[i * 3 + 1] + r22*T[i * 3 + 2] + t2; query[2] = sz;

			// search nearest neighbor
			M_tree->n_nearest(query, 1, neighbor);

			// model point
			double dx = M_tree->the_data[neighbor[0].idx][0];
			double dy = M_tree->the_data[neighbor[0].idx][1];
			double dz = M_tree->the_data[neighbor[0].idx][2];

			// model point normal
			double nx = M_normal[neighbor[0].idx * 3 + 0];
			double ny = M_normal[neighbor[0].idx * 3 + 1];
			double nz = M_normal[neighbor[0].idx * 3 + 2];

			// check if it is an inlier
			if ((sx - dx)*nx + (sy - dy)*ny + (sz - dz)*nz < indist)
			{
				fit_set.push_back(neighbor[0].idx);
			}
		}
	}
}

void IcpPointToPlane::computeNormal (const kdtree::KDTreeResultVector &neighbors,double *M_normal,const double flatness) {
  
  // dimensionality 2
  if (dim==2) {
    
    // extract neighbors
    Matrix P(neighbors.size(),2);
    Matrix mu(1,2);
    for (int i=0; i<neighbors.size(); i++) {
      double x = M_tree->the_data[neighbors[i].idx][0];
      double y = M_tree->the_data[neighbors[i].idx][1];
      P.val[i][0] = x;
      P.val[i][1] = y;
      mu.val[0][0] += x;
      mu.val[0][1] += y;
    }

    // zero mean
    mu       = mu/(double)neighbors.size();
    Matrix Q = P - Matrix::ones(neighbors.size(),1)*mu;

    // principal component analysis
    Matrix H = ~Q*Q;
    Matrix U,W,V;
    H.svd(U,W,V);

    // normal
    M_normal[0] = U.val[0][1];
    M_normal[1] = U.val[1][1];
  
  // dimensionality 3
  } else {
    
    // extract neighbors
    Matrix P(neighbors.size(),3);
    Matrix mu(1,3);
    for (int i=0; i<neighbors.size(); i++) {
      double x = M_tree->the_data[neighbors[i].idx][0];
      double y = M_tree->the_data[neighbors[i].idx][1];
      double z = M_tree->the_data[neighbors[i].idx][2];
      P.val[i][0] = x;
      P.val[i][1] = y;
      P.val[i][2] = z;
      mu.val[0][0] += x;
      mu.val[0][1] += y;
      mu.val[0][2] += z;
    }

    // zero mean
    mu       = mu/(double)neighbors.size();
    Matrix Q = P - Matrix::ones(neighbors.size(),1)*mu;

    // principal component analysis
    Matrix H = ~Q*Q;
    Matrix U,W,V;
    H.svd(U,W,V);

    // normal
    M_normal[0] = U.val[0][2];
    M_normal[1] = U.val[1][2];
    M_normal[2] = U.val[2][2];
  }
}

double* IcpPointToPlane::computeNormals (const int num_neighbors,const double flatness) {
  double *M_normal = (double*)malloc(M_tree->N*dim*sizeof(double));
  kdtree::KDTreeResultVector neighbors;
  for (int i=0; i<M_tree->N; i++) {
    M_tree->n_nearest_around_point(i,0,num_neighbors,neighbors);
    if (dim==2) computeNormal(neighbors,M_normal+i*2,flatness);
    else        computeNormal(neighbors,M_normal+i*3,flatness);
  }
  return M_normal;
}
