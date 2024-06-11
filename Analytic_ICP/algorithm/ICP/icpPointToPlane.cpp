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
double IcpPointToPlane::fitStep(double *T, const int T_num, MatrixXd &R, VectorXd &t, const vector<int> active, vector<int> &fit_set)
{
	int active_num = active.size();
	// kd tree query + result
	std::vector<float>         query(dim);
	kdtree::KDTreeResultVector result;

	// init matrix for point correspondences
	MatrixXd p_m(active_num, dim); // model
	MatrixXd p_t(active_num, dim); // template
	fit_set.clear();
	// dimensionality 2
	if (dim == 2) {

		// extract matrix and translation vector
		double r00 = R(0, 0); double r01 = R(0, 1);
		double r10 = R(1, 0); double r11 = R(1, 1);
		double t0 = t(0); double t1 = t(1);

		// init A and b
		MatrixXd A(active_num, 3);
		MatrixXd b(active_num, 1);

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
			A(i, 0) = ny*sx - nx*sy;
			A(i, 1) = nx;
			A(i, 2) = ny;
			b(i, 0) = nx*dx + ny*dy - nx*sx - ny*sy;
		}

		// linear least square matrices
		MatrixXd A_ = A.transpose()*A;
		MatrixXd b_ = A.transpose()*b;

		// solve linear system
		if (solve(A_, b_)) {

			// rotation matrix
			MatrixXd R_ = MatrixXd::Identity(2, 2);
			R_(0, 1) = -b_(0, 0);
			R_(1, 0) = +b_(0, 0);

			// orthonormalized rotation matrix
			JacobiSVD<MatrixXd> svd(R_, ComputeThinU | ComputeThinV);
			MatrixXd U, W, V;
			U = svd.matrixU();

			V = svd.matrixV();

			W = svd.singularValues();

			R_ = U*V.transpose();

			// translation vector
			VectorXd t_(2);
			t_(0) = b_(1, 0);
			t_(1) = b_(2, 0);

			// compose: R|t = R_|t_ * R|t
			R = R_*R;
			t = R_*t + t_;
			return max((R_ - MatrixXd::Identity(2,2)).norm(), t_.norm());
		}

		// dimensionality 3
	}
	else {

		// extract matrix and translation vector
		double r00 = R(0, 0); double r01 = R(0, 1); double r02 = R(0, 2);
		double r10 = R(1, 0); double r11 = R(1, 1); double r12 = R(1, 2);
		double r20 = R(2, 0); double r21 = R(2, 1); double r22 = R(2, 2);
		double t0 = t(0); double t1 = t(1); double t2 = t(2);

		// init A and b
		MatrixXd A(active_num, 6);
		MatrixXd b(active_num, 1);

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
			A(i, 0) = nz*sy - ny*sz;
			A(i, 1) = nx*sz - nz*sx;
			A(i, 2) = ny*sx - nx*sy;
			A(i, 3) = nx;
			A(i, 4) = ny;
			A(i, 5) = nz;
			b(i, 0) = nx*dx + ny*dy + nz*dz - nx*sx - ny*sy - nz*sz;
		}

		// linear least square matrices
		MatrixXd A_ = A.transpose()*A;
		MatrixXd b_ = A.transpose()*b;

		// solve linear system
		if (solve(A_, b_)) {

			// rotation matrix
			MatrixXd R_ = MatrixXd::Identity(3, 3);
			R_(0, 1) = -b_(2, 0);
			R_(1, 0) = +b_(2, 0);
			R_(0, 2) = +b_(1, 0);
			R_(2, 0) = -b_(1, 0);
			R_(1, 2) = -b_(0, 0);
			R_(2, 1) = +b_(0, 0);

			// orthonormalized rotation matrix
			MatrixXd U, W, V;
			JacobiSVD<MatrixXd> svd(R_, ComputeThinU | ComputeThinV);
			U = svd.matrixU();

			V = svd.matrixV();

			W = svd.singularValues();
			R_ = U*V.transpose();

			// translation vector
			VectorXd t_(3);
			t_(0) = b_(3, 0);
			t_(1) = b_(4, 0);
			t_(2) = b_(5, 0);

			// compose: R|t = R_|t_ * R|t
			R = R_*R;
			t = R_*t + t_;
			return max((R_ - MatrixXd::Identity(3,3)).norm(), t_.norm());
		}
	}

	// failure
	return 0;
}

void IcpPointToPlane::getInliers(vector<int> &fit_set, double *T, const int T_num
	, const MatrixXd &R, const VectorXd &t, const double indist) {

	// init inlier vector + query point + query result
	vector<int>            inliers;
	std::vector<float>         query(dim);
	kdtree::KDTreeResultVector neighbor;
	fit_set.clear();

	// dimensionality 2
	if (dim == 2) {

		// extract matrix and translation vector
		double r00 = R(0, 0); double r01 = R(0, 1);
		double r10 = R(1, 0); double r11 = R(1, 1);
		double t0 = t(0); double t1 = t(1);

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
		double r00 = R(0, 0); double r01 = R(0, 1); double r02 = R(0, 2);
		double r10 = R(1, 0); double r11 = R(1, 1); double r12 = R(1, 2);
		double r20 = R(2, 0); double r21 = R(2, 1); double r22 = R(2, 2);
		double t0 = t(0); double t1 = t(1); double t2 = t(2);

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
	MatrixXd P(neighbors.size(), 2);
	VectorXd mu(2);
    for (int i=0; i<neighbors.size(); i++) {
      double x = M_tree->the_data[neighbors[i].idx][0];
      double y = M_tree->the_data[neighbors[i].idx][1];
	  P(i, 0) = x;
	  P(i, 1) = y;
      mu(0) += x;
      mu(1) += y;
    }

    // zero mean
    mu       = mu/(double)neighbors.size();
	MatrixXd Q = P - MatrixXd::Ones(neighbors.size(), 1)*mu;

    // principal component analysis
	MatrixXd H = Q.transpose()*Q;
	MatrixXd U, W, V;
	JacobiSVD<MatrixXd> svd(H, ComputeThinU | ComputeThinV);
	U = svd.matrixU();

	V = svd.matrixV();

	W = svd.singularValues();

    // normal
	M_normal[0] = U(0, 1);
	M_normal[1] = U(1, 1);
  
  // dimensionality 3
  } else {
    
    // extract neighbors
	MatrixXd P(neighbors.size(), 3);
	VectorXd mu(3);
    for (int i=0; i<neighbors.size(); i++) {
      double x = M_tree->the_data[neighbors[i].idx][0];
      double y = M_tree->the_data[neighbors[i].idx][1];
      double z = M_tree->the_data[neighbors[i].idx][2];
	  P(i, 0) = x;
	  P(i, 1) = y;
	  P(i, 2) = z;
      mu(0) += x;
      mu(1) += y;
      mu(2) += z;
    }

    // zero mean
    mu       = mu/(double)neighbors.size();
	MatrixXd Q = P - MatrixXd::Ones(neighbors.size(), 1)*mu;

    // principal component analysis
	MatrixXd H = Q.transpose()*Q;
	MatrixXd U, W, V;
	JacobiSVD<MatrixXd> svd(H, ComputeThinU | ComputeThinV);
	U = svd.matrixU();

	V = svd.matrixV();

	W = svd.singularValues();

    // normal
	M_normal[0] = U(0, 2);
	M_normal[1] = U(1, 2);
	M_normal[2] = U(2, 2);
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


bool IcpPointToPlane::solve(MatrixXd &A, MatrixXd &B, double eps) 
{
	if (A.rows() != A.cols() || A.rows() != B.rows() || A.rows() < 1 || B.cols() < 1) {
		cerr << "ERROR: Trying to eliminate matrices of size (" << A.rows() << "x" << A.cols() <<
			") and (" << B.rows() << "x" << B.cols() << ")" << endl;
		exit(0);
	}
	int m = B.rows();
	int n = B.cols();
	// index vectors for bookkeeping on the pivoting
	int* indxc = new int[m];
	int* indxr = new int[m];
	int* ipiv = new int[m];

	// loop variables
	int i, icol, irow, j, k, l, ll;
	double big, dum, pivinv, temp;

	// initialize pivots to zero
	for (j = 0; j < m; j++)
		ipiv[j] = 0;

	// main loop over the columns to be reduced
	for (i = 0; i < m; i++) {

		big = 0.0;

		// search for a pivot element
		for (j = 0; j < m; j++)
			if (ipiv[j] != 1)
				for (k = 0; k < m; k++)
					if (ipiv[k] == 0)
						if (fabs(A(j, k)) >= big)
						{
							big = fabs(A(j, k));
							irow = j;
							icol = k;
						}
		++(ipiv[icol]);

		// We now have the pivot element, so we interchange rows, if needed, to put the pivot
		// element on the diagonal. The columns are not physically interchanged, only relabeled.
		if (irow != icol) {
			for (l = 0; l < m; l++)
			{
				SWAP(A(irow, l), A(icol, l))
			}
			for (l = 0; l < n; l++) 
			{ 
				SWAP(B(irow,l), B(icol,l))
			}
		}

		indxr[i] = irow; // We are now ready to divide the pivot row by the
		indxc[i] = icol; // pivot element, located at irow and icol.

		// check for singularity
		if (fabs(A(icol, icol)) < eps) {
			delete[] indxc;
			delete[] indxr;
			delete[] ipiv;
			return false;
		}

		pivinv = 1.0 / A(icol,icol);
		A(icol,icol) = 1.0;
		for (l = 0; l < m; l++) A(icol,l) *= pivinv;
		for (l = 0; l < n; l++) B(icol,l) *= pivinv;

		// Next, we reduce the rows except for the pivot one
		for (ll = 0; ll < m; ll++)
			if (ll != icol) {
				dum = A(ll,icol);
				A(ll,icol) = 0.0;
				for (l = 0; l < m; l++)
					A(ll,l) -= A(icol,l) * dum;
				for (l = 0; l < n; l++)
					B(ll,l) -= B(icol,l) * dum;
			}
	}

	// This is the end of the main loop over columns of the reduction. It only remains to unscramble
	// the solution in view of the column interchanges. We do this by interchanging pairs of
	// columns in the reverse order that the permutation was built up.
	for (l = m - 1; l >= 0; l--) {
		if (indxr[l] != indxc[l])
			for (k = 0; k < m; k++)
				SWAP(A(k,indxr[l]), A(k,indxc[l]))
	}

	// success
	delete[] indxc;
	delete[] indxr;
	delete[] ipiv;
	return true;
}