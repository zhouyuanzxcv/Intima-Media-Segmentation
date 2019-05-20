#define _USE_MATH_DEFINES

#include <string.h>
#include <math.h>
#include <vector>
#include "assert.h"
#include "config.h"


#if defined(_USE_SPM4D_HASH)
#include "spm4d_hash.h"
typedef spm4d_hash<double> spm4d_type;
#elif defined(_USE_SPM4D_TREE)
#include "spm4d_tree.h"
typedef spm4d_tree<double> spm4d_type;
#elif defined(_USE_SPM4D_TRIE)
#include "spm4d_trie.h"
typedef spm4d_trie<double> spm4d_type;
#endif /*defined(_USE_SPM4D_HASH)*/


using namespace std;

#define ACCESS_SEGMENT(a,j) a+m*(j)
#define f_seg_(j) ACCESS_SEGMENT(f,j)
#define grad_dir_seg_(j) ACCESS_SEGMENT(grad_dir,j)

#define calc_y_coord(rho,theta,x) ((rho)/sin(theta)-((x)+1)*tan(M_PI/2-(theta)))
#define left_end(rho,theta) calc_y_coord(rho,theta,-1)
#define right_end(rho,theta,width) calc_y_coord(rho,theta,width-1)

struct dp_cell {
  double theta1,theta2,rho1,rho2;
  double cost;
  int prev;
};

extern void cdld(double *f, double *grad_dir, int m, int n, double epsilon, 
    int rho_num, int theta_num, int dir_tol, double angle_max, double dis_min, double dis_max, 
    spm4d_type::cell_list &oA, double *orhos, double *othetas);

void get_line(double rho, double theta, double *y, int n);

void dld_dp(double *f, double *grad_dir, int m, int n, int seg_num, double epsilon, 
    double t1, double kappa, int rho_num, int theta_num, int dir_tol, double w1, double w2, 
    double dis_min, double dis_max, double *y1, double *y2) {
  //double epsilon = 0.2; // threshold for categorizing pixel to edge point
  //double t1 = 0.3; // threshold for removing small elements in the 4-D matrix of CDLD
  //double kappa = 0.5; // parameter before continuous constraint
  
  int seg_len = ceil((double)n/seg_num);
  double *rhos = new double[rho_num+1]; 
  double *thetas = new double[theta_num];

  typedef vector<dp_cell> dp_cell_vec;
  dp_cell_vec* dp_mat = new dp_cell_vec[seg_num]; // dynamic programming cost map

  double min_value,max_value;
  int start_col,end_col,width;
  for (int i = 0; i < seg_num; i++) {
    start_col = i*seg_len;
    end_col = GW_MIN((i+1)*seg_len,n);
    width = end_col-start_col;

    double* seg_img = f_seg_(start_col);
    double* seg_dir = grad_dir_seg_(start_col);

    // normalize image segment
    min_value = seg_img[0];
    max_value = seg_img[0];
    for (int j = 0; j < m*width; j++) {
      if (seg_img[j] < min_value) min_value = seg_img[j];
      if (seg_img[j] > max_value) max_value = seg_img[j];
    }
    double diff = max_value-min_value;
    for (int j = 0; j < m*width; j++) {
      seg_img[j] = (seg_img[j]-min_value)/diff;
    }

    // constrained dual line detection
    spm4d_type::cell_list A;
    cdld(seg_img, seg_dir, m, end_col-start_col, epsilon, rho_num, theta_num, dir_tol, 
        w1, dis_min, dis_max, A, rhos, thetas);
    
    // remove cells with small value and initialize DP cost map
    max_value = A.begin()->value;
    spm4d_type::cell_list::const_iterator pos;
    for (pos = A.begin(); pos != A.end(); pos++) {
      if (pos->value > max_value) max_value = pos->value;      
    }
    dp_cell dc;
    for (pos = A.begin(); pos != A.end(); pos++) {
      if (pos->value > t1*max_value) {
        dc.theta1 = thetas[pos->d1];
        dc.theta2 = thetas[pos->d2];
        dc.rho1 = rhos[pos->d3];
        dc.rho2 = rhos[pos->d4];
        dc.cost = -(pos->value);
        dc.prev = -1;
        dp_mat[i].push_back(dc);
      }
    }
    // calculate cost map
    double cost,new_cost;
    dp_cell *plc,*prc;
    double smooth1,smooth2;
    if (i > 0) {
      for (int j = 0; j < dp_mat[i].size(); j++) {
        cost = GW_INFINITE;
        prc = &dp_mat[i][j];

        for (int k = 0; k < dp_mat[i-1].size(); k++) {
          plc = &dp_mat[i-1][k];          

          if ( abs(plc->theta1 - prc->theta1) <= w2 && abs(plc->theta2 - prc->theta2) <= w2 ) {
            smooth1 = abs(right_end(plc->rho1,plc->theta1,seg_len)-left_end(prc->rho1,prc->theta1));
            smooth2 = abs(right_end(plc->rho2,plc->theta2,seg_len)-left_end(prc->rho2,prc->theta2));
            new_cost = (plc->cost)+kappa*(smooth1+smooth2);
            if (new_cost < cost) {
              cost = new_cost;
              prc->prev = k;
            }
          }
        }
        prc->cost += cost;
      }
    }
  }

  // back tracing
  assert(dp_mat[seg_num-1].size() > 0);

  int ka = 0;  
  min_value = dp_mat[seg_num-1][ka].cost;
  for (int i = 0; i < dp_mat[seg_num-1].size(); i++)
    if (dp_mat[seg_num-1][i].cost < min_value) {
      min_value = dp_mat[seg_num-1][i].cost;
      ka = i;
    }

  dp_cell *pcc;
  for (int i = seg_num-1; i >= 0; i--) {
    pcc = &dp_mat[i][ka];
    start_col = i*seg_len;
    end_col = GW_MIN((i+1)*seg_len,n);
    width = end_col-start_col;

    get_line(pcc->rho1,pcc->theta1,y1+start_col,width);
    get_line(pcc->rho2,pcc->theta2,y2+start_col,width);

    if (i > 0) ka = pcc->prev;
  }

  delete [] rhos;
  delete [] thetas;
  delete [] dp_mat;
}

void get_line(double rho, double theta, double *y, int n) {
  for (int i = 0; i < n; i++)
    y[i] = calc_y_coord(rho,theta,i);
}
