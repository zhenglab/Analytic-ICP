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

#ifndef ICP_POINT_TO_PLANE_H
#define ICP_POINT_TO_PLANE_H

#include "icp.h"

#define SWAP(a,b) {temp=a;a=b;b=temp;}

class IcpPointToPlane : public Icp {

public:
  
  IcpPointToPlane (double *M,const int M_num,const int dim,const int num_neighbors=10,const double flatness=5.0) : Icp(M,M_num,dim) {
    M_normal = computeNormals(num_neighbors,flatness);
  }

  IcpPointToPlane(Vector2d *&ps, int &p_num, const int num_neighbors = 10, const double flatness = 5.0) : Icp(ps, p_num) {
	  M_normal = computeNormals(num_neighbors, flatness);
  }

  IcpPointToPlane(Vector3d *&ps, int &p_num, const int num_neighbors = 10, const double flatness = 5.0) : Icp(ps, p_num) {
	  M_normal = computeNormals(num_neighbors, flatness);
  }

  virtual ~IcpPointToPlane () {
    delete M_normal;
  }

private:

	double fitStep(double *T, const int T_num, MatrixXd &R, VectorXd &t, const vector<int> active, vector<int> &fit_set);
	void getInliers(vector<int> &fit_set, double *T, const int T_num, const MatrixXd &R, const VectorXd &t, const double indist);
  
  // utility functions to compute normals from the model tree
  void computeNormal (const kdtree::KDTreeResultVector &neighbors,double *M_normal,const double flatness);
  double* computeNormals (const int num_neighbors,const double flatness);
  bool solve(MatrixXd &A, MatrixXd &B, double eps = 1e-20);
  // normals of model points
  double *M_normal;
};

#endif // ICP_POINT_TO_PLANE_H
