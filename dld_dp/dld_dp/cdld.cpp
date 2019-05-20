#define _USE_MATH_DEFINES

#include <string.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <list>
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

#define ACCESS_2D_ARRAY(a,i,j) a[(i)+m*(j)]
#define f_(i,j) ACCESS_2D_ARRAY(f,i,j)
#define grad_dir_(i,j) ACCESS_2D_ARRAY(grad_dir,i,j)

//#define mod(a,b) ((a) < 0 ? (fmod(double(a),(b))+(b)) : (fmod(double(a),(b))))
#define calc_dist(x,y,theta) (((x)+1)*cos(theta)+((y)+1)*sin(theta)) // matlab uses one based indexing


void cdld(double *f, double *grad_dir, int m, int n, double epsilon, int rho_num, int theta_num, int dir_tol, double angle_max, double dis_min, double dis_max, spm4d_type::cell_list &oA, double *orhos, double *othetas) {
  vector<int> ys,xs; 
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if (f_(i,j) > epsilon) {
        ys.push_back(i);
        xs.push_back(j);
      }

  double rho_min = 0;
  //double rho_max = sqrt(double(m*m+n*n));
  double rho_max = m;
  double theta_min = 0;
  double theta_max = M_PI;
  double rho_delta = (rho_max-rho_min)/rho_num;
  double theta_delta = (theta_max-theta_min)/theta_num;

  int M = rho_num+1;
  int N = theta_num;
  double *rhos = new double[M];
  double *thetas = new double[N];
  for (int i = 0; i < M; i++) {
    rhos[i] = rho_min + rho_delta*i;
  }
  for (int i = 0; i < N; i++) { 
    thetas[i] = theta_min + theta_delta*i;
  }

  spm4d_type A(N,N,M,M);  
  int angle_ind_max = GW_ROUND(angle_max/theta_delta);
  int dis_ind_min = GW_ROUND(dis_min/rho_delta);
  int dis_ind_max = GW_ROUND(dis_max/rho_delta);

//#if defined(_USE_SPM4D_HASH)
//  // set number of bucket
//  int bucket_num = 2*(xs.size())*(dir_tol*2+1)*(2*angle_ind_max+1)*(dis_ind_max-dis_ind_min+1);
//  // next bucket count obtained from the list of primes in the sgi hash_map  
//  bucket_num = stlp_std::priv::_Stl_prime_type::_S_next_size(bucket_num*0.12); /* 0.12 is an empirical value used to determine the number of buckets */
//  A.mat.resize(bucket_num); 
//#if defined(_DEBUG)
//  cout << "number of buckets before addition: " << bucket_num << endl;
//#endif
//#endif

  for (int i = 0; i < xs.size(); i++) {
    int y = ys[i];
    int x = xs[i];
    double dir = grad_dir_(y,x);
    
    if (dir < 0) continue;

    int thetam_ind = GW_ROUND(dir/theta_delta);

    double theta1,theta2,rho1,rho2;
    int theta1_ind,theta2_ind,rho1_ind,rho2_ind;
    for (int tol_ind = GW_MAX(thetam_ind-dir_tol,0); tol_ind <= GW_MIN(thetam_ind+dir_tol,N-1); tol_ind++) {      
      // row
      theta1_ind = tol_ind;
      theta1 = thetas[theta1_ind];
      rho1 = calc_dist(x,y,theta1);
      rho1_ind = GW_ROUND(rho1/rho_delta);     
      if (rho1_ind >= 0 && rho1_ind < M) {
        for (theta2_ind = GW_MAX(theta1_ind-angle_ind_max,0); theta2_ind <= GW_MIN(theta1_ind+angle_ind_max,N-1); theta2_ind++) 
          for (rho2_ind = GW_MAX(rho1_ind+dis_ind_min,0); rho2_ind <= GW_MIN(rho1_ind+dis_ind_max,M-1); rho2_ind++)
            A.increment_by(theta1_ind,theta2_ind,rho1_ind,rho2_ind,f_(y,x));        
      }
      // column
      theta2_ind = tol_ind;
      theta2 = thetas[theta2_ind];
      rho2 = calc_dist(x,y,theta2);
      rho2_ind = GW_ROUND(rho2/rho_delta);
      if (rho2_ind >= 0 && rho2_ind < M) {
        for (theta1_ind = GW_MAX(theta2_ind-angle_ind_max,0); theta1_ind <= GW_MIN(theta2_ind+angle_ind_max,N-1); theta1_ind++) 
          for (rho1_ind = GW_MAX(rho2_ind-dis_ind_max,0); rho1_ind <= GW_MIN(rho2_ind-dis_ind_min,M-1); rho1_ind++)
            A.increment_by(theta1_ind,theta2_ind,rho1_ind,rho2_ind,f_(y,x));        
      }
    }
  }

  // output
  A.convert_to_list(oA);

//#if defined(_USE_SPM4D_HASH) && defined(_DEBUG)
//  cout << "number of elements inserted: " << oA.size() << endl;
//  cout << "number of buckets after addition: " << A.mat.bucket_count() << endl << endl;
//#endif

  for (int i = 0; i < M; i++) orhos[i] = rhos[i];
  for (int i = 0; i < N; i++) othetas[i] = thetas[i];    

  delete [] rhos;
  delete [] thetas;  
}
