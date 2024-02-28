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

#ifndef ICP_POINT_TO_POINT_H
#define ICP_POINT_TO_POINT_H

#include "icp.h"

class IcpPointToPoint : public Icp {

public:

  IcpPointToPoint (double *M,const int M_num,const int dim) : Icp(M,M_num,dim) {}
  IcpPointToPoint(CPoint2 *&ps, int &p_num) : Icp(ps, p_num) {}
  IcpPointToPoint(CPoint2 *&ps, int &p_num, int iterCount) : Icp(ps, p_num, iterCount) {}

  IcpPointToPoint(CPoint3 *&ps, int &p_num) : Icp(ps, p_num) {}
  IcpPointToPoint(CPoint3 *&ps, int &p_num, int iterCount) : Icp(ps, p_num, iterCount) {}
  virtual ~IcpPointToPoint () {}

private:

  double fitStep (double *T,const int T_num,Matrix &R,Matrix &t, const vector<int> active, vector<int> &fit_set);
  void getInliers (vector<int> &fit_set, double *T,const int T_num,const Matrix &R,const Matrix &t,const double indist);
};

#endif // ICP_POINT_TO_POINT_H
