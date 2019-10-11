#ifndef BA_HPP
#define BA_HPP

#include <vector>
#include <stdio.h>
#include <Eigen/Eigen>
#include <Eigen/Dense>

#include <math.h>

#define NOT_FOUND -1
typedef std::vector<std::vector<int> > Table;
struct Point3f{
  float x, y, z;
};
struct Camera {
  float r1, r2, r3, t1, t2, t3, f, k1,k2;
  public:
  float operator[](int i){
    if(i == 0)
      return r1;
    if(i == 1) 
      return r2;
    if(i == 2) 
      return r3;
    if(i == 3) 
      return t1;
    if(i == 4) 
      return t2;
    if(i == 5) 
      return t3;
    if(i == 6) 
      return f;
    if(i == 7) 
      return k1;
    if(i == 8) 
      return k2;
  }
};
struct Observation{
  int cid;
  int pid;
  float x;
  float y;
};
struct BAProblem{
  int camera_number;
  int point_number;
  int observation_number;
  //if x* is a good approximation tau can be set as a small
  //float like 1.0e-6 otherwise set tau as 1e-3 or even 1.0

  double tau;
  std::vector<Point3f> point_data;
  std::vector<Camera> camera_data;
  std::vector<Observation> observation_data;
  //delta a delta b 
  //delta = [delta_a; delta_b], delta should be calculated 
  //using (A+lambda*I)*delta = g;
  //p stores the optimized parameters
  // p_new = p_old + delta
  Eigen::VectorXd delta_a;  
  Eigen::VectorXd delta_b;  
  Eigen::VectorXd delta;
  int k;
  int v;
  double tau;
  

};

#endif
