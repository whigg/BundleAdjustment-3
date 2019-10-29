#ifndef BA_HPP
#define BA_HPP

#include <vector>
#include <stdio.h>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <iostream>
#include <math.h>

#define NOT_FOUND -1
typedef std::vector<std::vector<int> > Table;
typedef std::vector<int> vint;

void rodrigues(Eigen::Vector3d const rvec, Eigen::Matrix3d &rotm) {
  double theta = rvec.norm();
  Eigen::Vector3d r = rvec;
  r = r / theta;
  Eigen::Matrix<double,3,3> I = Eigen::Matrix<double, 3,3>::Identity();
  Eigen::Matrix<double,3,3> rx;
  rx << 0.0, -r(2), r(1), r(2), 0.0, -r(0), -r(1), r(0), 0.0;
  rotm = cos(theta) * I + (1-cos(theta))*r*r.transpose() + sin(theta)* rx;
}

struct Point3f{
  float x, y, z;
};
struct Camera {
  float r1, r2, r3, t1, t2, t3, f, k1,k2;
  public:
  float operator[](int i){
    if(i>7 || i <0){
      printf("index exceeded\n");
      exit(-1);
    }
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
  private:
    int camera_number;
    int point_number;
    int observation_number;
    //if x* is a good approximation tau can be set as a small
    //float like 1.0e-6 otherwise set tau as 1e-3 or even 1.0

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
    Eigen::VectorXd g;
    Eigen::VectorXd p;//camera parameter and point coord
    std::vector<int> base;//index of observation_data
    Table captt; //camera point table. captt[cid] return corresponding pid.

    int k;
    int v;
    double tau;
  public:
    //reprojection error of point pid projected to cid
    double repj_err(int const pid, int const cid)const {
	int obs_i = this->at(pid, cid);
	if(obs_i == NOT_FOUND) {
	  return 0.0;
	}
	float obsx, obsy;
	obsx = observation_data[obs_i].x;
	obsy = observation_data[obs_i].y;
	Camera cam = camera_data[cid];
	Point3f pt = point_data[pid];
	Eigen::Vector3d rvec;
	rvec(0,0) = cam.r1;
	rvec(1,0) = cam.r2;
	rvec(2,0) = cam.r3;
	float f = cam.f;
	float k1 = cam.k1;
	float k2 = cam.k2;
	Eigen::Vector3d t;
	t(0) = cam.t1;
	t(1) = cam.t2;
	t(2) = cam.t3;
	Eigen::Matrix3d rotm;
	rodrigues(rvec, rotm);
	Eigen::Vector3d Xw ;
	Xw << pt.x, pt.y, pt.z;
	Eigen::Vector3d Xc = rotm * Xw + t;
	Xc = -Xc/Xc(2);
	double l1 = Xc.norm();
	double l2 = l1*l1;
	double rho = 1.0 + k1 * l1*l1 + k2 * l2*l2;
	printf("rho = %lf\n", rho);
	Eigen::Vector3d Xp = f*rho*Xc;
	float d1 = fabs(Xp(0,0) - obsx);
	float d2 = fabs(Xp(1,0) - obsy);
	return sqrt(d1*d1+d2*d2);	  
    }
  
    int at(int const pid, int const cid) const {
      if(pid >= point_number || cid >= camera_number) {
      	printf("point index or camera index exceeded\n");
        exit(-1);
      }
      int base = this->base[pid];
      int limit = this->base[pid+1];
      for(int i=base; i<limit; i++) {
        Observation obs;
        obs = observation_data[i];
        if(pid == obs.pid && cid == obs.cid) {
          return i;
        }
      }
      return NOT_FOUND;
    }
    void at(int const pid, vint & vcid )const {
      if(pid >= point_number) {
        printf("pid index exceeded\n");
        exit(-1);
      }
      int base = this->base[pid];
      int limit = this->base[pid+1];
      vcid.clear();
      for(int i=base; i<limit; i++) {
        Observation obs;
        obs = observation_data[i];
        vcid.push_back(obs.cid);
      }
    }
    BAProblem(char const *file_name) {
      tau = 1.0e-4;
      FILE *fp = fopen(file_name, "r");
      if(!fp) {perror("error opening file\n"); exit(-1);}
      fscanf(fp, "%d", &camera_number);
      fscanf(fp, "%d", &point_number);
      fscanf(fp, "%d", &observation_number);
      vint blocks;
      blocks.clear();
      this->base.clear();
      int current_pid = 0;
      int current_pid_num = 0;
      for(int i=0;i<camera_number; i++) {
        vint empty_vec;
        empty_vec.clear();
        captt.push_back(empty_vec);
      }
      for(int i=0; i<observation_number; i++) {
        Observation obs;
        fscanf(fp, "%d", &(obs.cid));
        fscanf(fp, "%d", &(obs.pid));
        fscanf(fp, "%f", &(obs.x));
        fscanf(fp, "%f", &(obs.y));
        captt[obs.cid].push_back(obs.pid);
        if(current_pid == obs.pid) {
          current_pid_num++;
        }else {
          blocks.push_back(current_pid_num);
          current_pid_num = 1;
	        current_pid = obs.pid;
        }
        observation_data.push_back(obs);
      }
      blocks.push_back(current_pid_num);
      base.push_back(0);
      for(unsigned i=0; i<blocks.size(); i++) {
        base.push_back(base[i]+blocks[i]);
      }
      camera_data.clear();
      p = Eigen::VectorXd(camera_number*6 + point_number*3);
      for(int i=0; i<camera_number; i++) {
        Camera cam;
        fscanf(fp, "%f", & (cam.r1));
        fscanf(fp, "%f", & (cam.r2));
        fscanf(fp, "%f", & (cam.r3));
        fscanf(fp, "%f", & (cam.t1));
        fscanf(fp, "%f", & (cam.t2));
        fscanf(fp, "%f", & (cam.t3));
        fscanf(fp, "%f", & (cam.f ));
        fscanf(fp, "%f", & (cam.k1));
        fscanf(fp, "%f", & (cam.k2));
        camera_data.push_back(cam);
        for(int j=0; j<6; j++) {
          p(i*6+j) = cam[j];
        }
      }
      point_data.clear();
      for(int i=0; i<point_number; i++) {
        Point3f pt;
        fscanf(fp, "%f", &(pt.x));
        fscanf(fp, "%f", &(pt.y));
        fscanf(fp, "%f", &(pt.z));
        point_data.push_back(pt);
        int b= camera_number*6;
        p(b + i*3+0) = pt.x; p(b+i*3+1) = pt.y; p(b+i*3+2) = pt.z;
      }
      fclose(fp);
    }
    void test() {
      vint a = this->captt[49];
      for(auto it=a.begin(); it!= a.end(); it++) {
        printf("%d ", *it);
      }
      printf("\n");
    }

};

#endif
