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

#ifndef ICP_H
#define ICP_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

#include "matrix.h"
#include "kdtree.h"

using namespace std;

class Icp
{

public:

  // constructor
  // input: M ....... pointer to first model point
  //        M_num ... number of model points
  //        dim   ... dimensionality of model points (2 or 3)
  Icp(double *M, const int M_num, const int dim);
  Icp(CPoint2 *&ps, int &p_num, int iterCount);
  Icp(CPoint2 *&ps, int &p_num);

  Icp(CPoint3 *&ps, int &p_num, int iterCount);
  Icp(CPoint3 *&ps, int &p_num);
  
  // deconstructor
  virtual ~Icp ();
  
  // set subsampling step size of coarse matching (1. stage)
  void setSubsamplingStep (int val) { sub_step  = val; }
  
  // set maximum number of iterations (1. stopping criterion)
  void setMaxIterations   (int val) { max_iter  = val; }
  
  // set minimum delta of rot/trans parameters (2. stopping criterion)
  void setMinDeltaParam   (double  val) { min_delta = val; }
  
  // fit template to model yielding R,t (M = R*T + t)
  // input:  T ....... pointer to first template point
  //         T_num ... number of template points
  //         R ....... initial rotation matrix
  //         t ....... initial translation vector
  //         indist .. inlier distance (if <=0: use all points)
  // output: R ....... final rotation matrix
  //         t ....... final translation vector
  void fit(double *T,const int T_num,Matrix &R,Matrix &t, vector<int> &fit_set, const double indist);

  void fit(CPoint2 *&ps, int &p_num, Matrix &R, Matrix &t, vector<int> &fit_set, const double indist);
  
  void fit(CPoint3 *&ps, int &p_num, Matrix &R, Matrix &t, vector<int> &fit_set, const double indist);

private:
  
  // iterative fitting
  void fitIterate(double *&T, const int &T_num,Matrix &R,Matrix &t, const vector<int> active, vector<int> &fit_set);
  
  // inherited classes need to overwrite these functions
  virtual double fitStep(double *T,const int T_num,Matrix &R,Matrix &t, const vector<int> active, vector<int> &fit_set) = 0;
  virtual void getInliers(vector<int> &fit_set, double *T, const int T_num, const Matrix &R, const Matrix &t, const double indist) = 0;
  
protected:
  
  // kd tree of model points
  kdtree::KDTree*     M_tree;
  kdtree::KDTreeArray M_data;
  
  int dim;       // dimensionality of model + template data (2 or 3)
  int sub_step;  // subsampling step size
  int max_iter;  // max number of iterations
  double  min_delta; // min parameter delta

  void pointset2T(double *T, CPoint2 *ps, const int p_num);
  void pointset3T(double *T, CPoint3 *ps, const int p_num);
};

#endif // ICP_H
