/*=================================================================
% linked dual line detection using Hough transform and dynamic programming
%
%   [y1,y2] = dld_dp_mex(f,grad_dir,seg_num,epsilon,t1,kappa,rho_num,theta_num,dir_tol,w1,w2,dis_min,dis_max);
%
%   Input:
%     f - edge map
%     grad_dir - edge direction
%     seg_num - number of image segments
%     epsilon - threshold for determining the edge points
%     t1 - threshold for removing small values in the 4-D sparse matrix
%     kappa - curve smoothness (sum of squared distances between adjacent endpoints)
%     rho_num - number of discrete values for rho (delta_rho = M/rho_num)
%     theta_num - number of discrete values for theta (delta_theta = pi/theta_num)
%     dir_tol - tolerance of edge direction (dir_tol = theta_epsilon/delta_theta)
%     w1 - maximal angle between the two line segments in an image segment
%     w2 - maximal angle between adjacent two line segments 
%     dis_min - minimal distance between the two line segments in an image segment
%     dis_max - maximal distance between the two line segments in an image segment
%   Output:
%     y1 - vector of y-coordinates of boundary points of upper curve
%     y2 - vector of y-coordinates of boundary points of lower curve
%
%   Reference: Yuan Zhou, Xinyao Cheng, Xiangyang Xu, Enmin Song. Dynamic Programming
%   in Parallel Boundary Detection with Application to Ultrasound Intima-media Segmentation.
%   Medical Image Analysis, 2013
%   
%   The code uses SGI STL which has a faster implementation of hash_map compared to the version
%   that comes with Visual Studio. The SGI STL for Visual Studio 2010 can be found at 
%   http://www.lenholgate.com/blog/2010/07/stlport-521-and-vs2010-and-x64.html.
%
%   Copyright (c) 2012 ZHOU Yuan
*=================================================================*/

#include "mex.h"

extern void dld_dp(double *f, double *grad_dir, int m, int n, int seg_num, double epsilon, 
    double t1, double kappa, int rho_num, int theta_num, int dir_tol, double w1, double w2, 
    double dis_min, double dis_max, double *y1, double *y2);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) 
{ 
  // variables
  double *f;
  double *grad_dir;
  const mwSize *idims;
  mwSize odims[] = {1,1};
  int m, n; // dimension of edge map
  int seg_num;
  double epsilon,t1,kappa;
  int rho_num,theta_num,dir_tol;  
  double w1,w2,dis_min,dis_max;
  
  // parse input
  idims = mxGetDimensions(prhs[0]);
  m = (int)idims[0]; n = (int)idims[1];
  f = mxGetPr(prhs[0]);
  grad_dir = mxGetPr(prhs[1]);
  
  seg_num = (int)mxGetScalar(prhs[2]);

  epsilon = mxGetScalar(prhs[3]);
  t1 = mxGetScalar(prhs[4]);
  kappa = mxGetScalar(prhs[5]);

  rho_num = (int)mxGetScalar(prhs[6]);
  theta_num = (int)mxGetScalar(prhs[7]);
  dir_tol = (int)mxGetScalar(prhs[8]);

  w1 = mxGetScalar(prhs[9]);
  w2 = mxGetScalar(prhs[10]);
  dis_min = mxGetScalar(prhs[11]);
  dis_max = mxGetScalar(prhs[12]);
  
  // create output
  odims[0] = n;
  plhs[0] = mxCreateNumericArray(2, odims, mxDOUBLE_CLASS, mxREAL);
  plhs[1] = mxCreateNumericArray(2, odims, mxDOUBLE_CLASS, mxREAL);
  double *y1 = mxGetPr(plhs[0]);
  double *y2 = mxGetPr(plhs[1]);

  // compute dual line detection with dynamic programming
  dld_dp(f, grad_dir, m, n, seg_num, epsilon, t1, kappa, rho_num, theta_num, dir_tol, 
      w1, w2, dis_min, dis_max, y1, y2);
}
