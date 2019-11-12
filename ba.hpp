#ifndef BA_HPP
#define BA_HPP

#include <vector>
#include <stdio.h>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <iostream>
#include <math.h>
#include <assert.h>
#include <time.h>

#define NOT_FOUND -1
#define Delta 0.0005

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
  Eigen::MatrixXd delta_a;  
  Eigen::MatrixXd delta_b;
  Eigen::MatrixXd delta;
  Eigen::MatrixXd ga;
  
  std::vector<int> base;//index of observation_data
  Table captt; //camera point table. captt[cid] return corresponding pid.

  std::vector<Eigen::Matrix<double,6,6> > U;
  std::vector<Eigen::Matrix<double,3,3> > V;
  std::vector<Eigen::Matrix<double,6,3> > W;
  std::vector<Eigen::Matrix<double,6,1> > EpsilonCam;//epsilon a
  std::vector<Eigen::Matrix<double,3,1> > EpsilonPt; //epsilon b
  Eigen::MatrixXd S;
  
  int k;
  int v;
  double tau;
  double u;
public:
  // calculate EpsilonCam
  // calculate U and V
  void calcU(void) {
    std::vector<int> pid_vec;
    EpsilonCam.clear();
    for(int cid=0; cid<camera_number; cid++) {
      Eigen::MatrixXd Us = Eigen::MatrixXd::Zero(6,6);
      Eigen::MatrixXd Ui = Eigen::MatrixXd::Zero(1,6);
      Eigen::MatrixXd epsa = Eigen::MatrixXd::Zero(6,1);
      pid_vec.clear();
      pid_vec = captt[cid];      
      for(auto it = pid_vec.begin(); it!=pid_vec.end(); it++){
	int pid = *it;
	Ui = this->Jacc(pid, cid);
	double e = this->repj_err(pid, cid);
	epsa += Ui.transpose() * e;
	Us += Ui.transpose() * Ui;
      }
      double m = maxd(Us);
      u = u < m? m:u;
      EpsilonCam.push_back(epsa);
      U.push_back(Us);
    }    
  }

  void calcV(void) {
    std::vector<int> cid_vec;
    EpsilonPt.clear();
    for(int pid = 0; pid < point_number; pid++) {
      Eigen::MatrixXd Vs = Eigen::MatrixXd::Zero(3,3);
      Eigen::MatrixXd Vi = Eigen::MatrixXd::Zero(1,3);
      Eigen::MatrixXd epsb=Eigen::MatrixXd::Zero(3,1);
      cid_vec.clear();
      this->at(pid, cid_vec);
      for(auto it = cid_vec.begin(); it != cid_vec.end(); it++) {
	int cid = *it;
	Vi = this->Jacp(pid, cid);
	Vs += Vi.transpose() * Vi;
	double e = this->repj_err(pid, cid);
	epsb += Vi.transpose() * e;
      }
      double m = maxd(Vs);
      u = u < m? m:u;
      EpsilonPt.push_back(epsb);
      V.push_back(Vs);
    }
    
  }
  void calcW(void) {//camera_number*point_number matrix
    W.clear();
    
    Eigen::MatrixXd Aij = Eigen::MatrixXd::Zero(1,6);
    Eigen::MatrixXd Bij = Eigen::MatrixXd::Zero(1,3);
    for(int cid = 0; cid<camera_number; cid++) {
      Eigen::MatrixXd wij = Eigen::MatrixXd::Zero(6,3);
      for(int pid = 0; pid<point_number; pid++) {
	if(this->at(pid,cid)==NOT_FOUND) {
	  W.push_back(wij);
	}else{
	  Aij = Jacc(pid, cid);
	  Bij = Jacp(pid, cid);
	  wij = Aij.transpose() * Bij;
	  W.push_back(wij);
	}
	
      }
    }
    
  }
  Eigen::MatrixXd getW(int const pid, int const cid) const {
    //Eigen::MatrixXd wij = Eigen::MatrixXd::Zero(6,3);
    return W[cid * point_number + pid];
  }

  //Calculate S
  void calcS(){
    //calcU();
    //calcV();
    //calcW();
    
    u = tau*u;
    
    S = Eigen::MatrixXd::Zero(6*camera_number, 6*camera_number);
    for(int j=0; j<camera_number; j++) {
      
      for(int k=0; k<camera_number; k++){
	S.block<6,6>(6*j, 6*k) = Eigen::MatrixXd::Zero(6,6);
	Eigen::Matrix3d Vs;
	for(int pid = 0; pid<point_number; pid++) {
	  if(at(pid,j) != NOT_FOUND && at(pid,k)!=NOT_FOUND) {
	    Vs = V[pid] + u * Eigen::MatrixXd::Identity(3,3);
	    S.block<6,6>(6*j,6*k) -= getW(pid,j)*Vs.inverse()*getW(pid,k).transpose();
     	  }
	}//endfor
	if(j == k) {
	  
	  S.block<6,6>(6*j,6*k) = (U[j]+u*Eigen::MatrixXd::Identity(6,6))+S.block<6,6>(6*j,6*k);
	}	
      }//endfor
      
    }//endfor
    
  }//endfunc
  
  //Calculate g, ga, gb
  void calcG(void) {
    ga = Eigen::MatrixXd(camera_number*6,1);
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(3,3);
    
    for(int cid = 0; cid < camera_number; cid++) {
      Eigen::VectorXd epsa = EpsilonCam[cid];
      Eigen::VectorXd s   = Eigen::MatrixXd::Zero(6,1);
      std::vector<int> pid_vec = captt[cid];
      for(auto it = pid_vec.begin(); it!= pid_vec.end(); it++) {
	int pid = *it;
	Eigen::MatrixXd wij = getW(pid, cid);
	Eigen::MatrixXd epsb = EpsilonPt[pid];
	s += wij * (V[pid]+u * I).inverse() * epsb;
      }
      ga.block<6,1>(cid*6,0) =  epsa - s;
    }
  }
  
  //reprojection error of point pid projected to cid
  double repj_err(int const pid, int const cid)const {
    int obs_i = this->at(pid, cid);
    if(obs_i == NOT_FOUND) {
      return 0.0;
    }
    float obsx, obsy;
    obsx = observation_data[obs_i].x;
    obsy = observation_data[obs_i].y;
    //printf("x,y = (%f,%f)\n",obsx, obsy);
    Camera cam = camera_data[cid];
    Point3f pt = point_data[pid];
    //printf("3dpt = %f, %f, %f\n",pt.x, pt.y, pt.z);
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
    // printf("camera parameter = %lf, %lf, %lf, %lf, %lf, %lf, %f, %f, %f\n",
    //        rvec(0,0), rvec(1,0), rvec(2,0), t(0), t(1), t(2),
    //        f, k1,k2);
	
    Eigen::Matrix3d rotm;
    rodrigues(rvec, rotm);
    Eigen::Vector3d Xw ;
    Xw << pt.x, pt.y, pt.z;
    Eigen::Vector3d Xc = rotm * Xw + t;
    Xc = -Xc/Xc(2);
    double l1 = Xc.norm();
    double l2 = l1*l1;
    double rho = 1.0 + k1 * l1*l1 + k2 * l2*l2;
    //printf("rho = %lf\n", rho);
    Eigen::Vector3d Xp = f*rho*Xc;
    double d1 = fabs(Xp(0,0) - obsx);
    double d2 = fabs(Xp(1,0) - obsy);
    return sqrt(d1*d1+d2*d2);	  
  }
  double repj_err(Camera const &cam, Point3f const &pt, Observation const & obs) const {
    Eigen::Vector3d rvec;
    rvec(0) = cam.r1;
    rvec(1) = cam.r2;
    rvec(2) = cam.r3;
    
    Eigen::Vector3d t;
    t(0) = cam.t1;
    t(1) = cam.t2;
    t(2) = cam.t3;
    
    double f = cam.f;
    double k1 = cam.k1;
    double k2 = cam.k2;
    Eigen::Matrix<double, 3,3> rotm;
    rodrigues(rvec, rotm);//rodrigues transformation
    Eigen::Matrix<double, 3,1> Xw;
    Xw(0,0) = pt.x;
    Xw(1,0) = pt.y;
    Xw(2,0) = pt.z;
    Eigen::Matrix<double,3,1> Xc = rotm * Xw + t;
    Xc = -Xc / Xc(2,0);
    double l1 = Xc.norm();
    double l2 = l1*l1;
    double rho = 1.0 + k1*l1*l1 + k2 * l2 * l2;
    Eigen::Matrix<double, 3,1> Xp = f * rho * Xc;
    double d1 = fabs(Xp(0,0)-obs.x);
    double d2 = fabs(Xp(1,0)-obs.y);
    return sqrt(d1*d1+d2*d2);
  }
  
  //camera jacobian matrix at(pid, cid), return 1*6
  Eigen::MatrixXd Jacc(int const pid, int const cid) const {
    if(pid >= point_number || cid >= camera_number) {
      printf("index exceeded.\n");
      exit(-1);
    }
    Eigen::MatrixXd Jcam = Eigen::MatrixXd::Zero(1,6);    
    int obs_i = this->at(pid, cid);
    if(obs_i == NOT_FOUND) {
      return Jcam;
    }
    Camera cam;
    Point3f pt;
    Observation obs;
    cam = camera_data[cid];
    pt = point_data[pid];    
    obs.pid = pid;
    obs.cid = cid;
    obs.x = observation_data[obs_i].x;
    obs.y = observation_data[obs_i].y;
    Camera cam2 = cam;
    Camera cam1 = cam;
    cam2.r1 += Delta;
    cam1.r1 -= Delta;
    double pj1 = repj_err(cam1, pt, obs);
    double pj2 = repj_err(cam2, pt, obs);
    Jcam(0,0) = (pj2 - pj1)/(Delta*2);
    cam2.r1 = cam.r1;
    cam1.r1 = cam.r1;

    cam2.r2 += Delta;
    cam1.r2 -= Delta;
    pj1 = repj_err(cam1, pt, obs);
    pj2 = repj_err(cam2, pt, obs);
    Jcam(0,1) = (pj2 - pj1)/(2*Delta);
    cam2.r2 = cam.r2;
    cam1.r2 = cam.r2;

    cam2.r3 += Delta;
    cam1.r3 -= Delta;

    pj1 = repj_err(cam1, pt, obs);
    pj2 = repj_err(cam2, pt, obs);
    Jcam(0,2) = (pj2 - pj1)/(2*Delta);
    cam2.r3 = cam.r3;
    cam1.r3 = cam.r3;

    //calculate Jt
    cam2.t1 += Delta;
    cam1.t1 -= Delta;

    pj1 = repj_err(cam1, pt, obs);
    pj2 = repj_err(cam2, pt, obs);
    Jcam(0,3) = (pj2 - pj1)/(2*Delta);
    cam2.t1 = cam.t1;
    cam1.t1 = cam.t1;

    cam2.t2 += Delta;
    cam1.t2 -= Delta;

    pj1 = repj_err(cam1, pt, obs);
    pj2 = repj_err(cam2, pt, obs);
    Jcam(0,4) = (pj2 - pj1)/(2*Delta);
    cam2.t2 = cam.t2;
    cam1.t2 = cam.t2;

    cam2.t3 += Delta;
    cam1.t3 -= Delta;
    pj1 = repj_err(cam1, pt, obs);
    pj2 = repj_err(cam2, pt, obs);
    Jcam(0,5) = (pj2 - pj1)/(2*Delta);
    cam2.t3 = cam.t3;
    cam1.t3 = cam.t3;
    return Jcam;
  }
  Eigen::MatrixXd Jacp(int const pid, int const cid) const {
    if(pid >= point_number || cid >= camera_number) {
      printf("index exceeded\n");
      exit(-1);
    }
    Camera cam;
    Point3f pt;
    Observation obs;
    Eigen::MatrixXd Jpt = Eigen::MatrixXd::Zero(1,3);
    int obs_i = this->at(pid, cid);
    if(obs_i == NOT_FOUND) {
      return Jpt;
    }
    obs.pid = pid;
    obs.cid = cid;
    obs.x = observation_data[obs_i].x;
    obs.y = observation_data[obs_i].y;
    pt = point_data[pid];
    cam = camera_data[cid];
    Point3f pt2 = pt;
    Point3f pt1 = pt;
    pt2.x += Delta;
    pt1.x -= Delta;
    double pj1 = repj_err(cam, pt1, obs);
    double pj2 = repj_err(cam, pt2, obs);
    Jpt(0,0) = (pj2 - pj1)/(2*Delta);
    pt2.x = pt.x;
    pt1.x = pt.x;

    pt2.y += Delta;
    pt1.y -= Delta;
    pj1 = repj_err(cam, pt1, obs);
    pj2 = repj_err(cam, pt2, obs);
    Jpt(0,1) = (pj2 - pj1)/(2*Delta);
    pt2.y = pt.y;
    pt1.y = pt.y;

    pt2.z += Delta;
    pt1.z -= Delta;
    pj1 = repj_err(cam, pt1, obs);
    pj2 = repj_err(cam, pt2, obs);
    Jpt(0,2) = (pj2 - pj1)/(2*Delta);
    pt2.z = pt.z;
    pt1.z = pt.z;
    return Jpt;    
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
  double maxd(Eigen::MatrixXd const diag) const {
    int m = diag.rows() <= diag.cols() ?diag.rows() : diag.cols();
    double maximum = diag(0,0);
    for(int i=0; i<m; i++){
      maximum = maximum > diag(i,i) ? maximum:diag(i,i);
    }
    return maximum;
  }
  
  BAProblem(char const *file_name) {
    tau = 1.0e-5;
    u = 0.0;
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
      
    }
    point_data.clear();
    for(int i=0; i<point_number; i++) {
      Point3f pt;
      fscanf(fp, "%f", &(pt.x));
      fscanf(fp, "%f", &(pt.y));
      fscanf(fp, "%f", &(pt.z));
      point_data.push_back(pt);
      //int b= camera_number*6;
      //p(b + i*3+0) = pt.x; p(b+i*3+1) = pt.y; p(b+i*3+2) = pt.z;
    }
    fclose(fp);
  }
  void solve_lin_eq(){
    //calcS();
    //calcG();//solve ga
    delta_a = S.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(ga); // solve delta a
    //solve delta_b
    delta_b = Eigen::MatrixXd::Zero(point_number*3,1);
    
    for(int pid = 0; pid < point_number; pid++) {
      Eigen::MatrixXd sum = Eigen::MatrixXd::Zero(3,1);
      Eigen::MatrixXd I = u * Eigen::MatrixXd::Identity(3,3);
      std::vector<int> cid_vec;
      this->at(pid, cid_vec);
      for(auto it=cid_vec.begin(); it!=cid_vec.end(); it++) {
	int cid = *it;
	Eigen::MatrixXd wij = getW(pid, cid);
	sum+= wij.transpose() * delta_a.block<6,1>(cid*6,0);
      }//endfor
      delta_b.block<3,1>(pid*3,0) = (V[pid]+I).inverse() * (EpsilonPt[pid] - sum);
    }
    
    //release();
    
  }
  
  void solve_lm(){
    //initialize step
    int k = 0;
    int max_iter = 10;
    delta_a = Eigen::MatrixXd::Zero(camera_number*6,1);
    Eigen::MatrixXd delta_new = delta_a;

    do{
      std::cout<<"iter "<<k<<std::endl;
      delta_new = delta_a;
      calcU();
      calcV();
      calcW();
      calcS();
      calcG();
      solve_lin_eq();
      update();
      release();
      k++;
      
      std::cout<<"norm = "<<(delta_new-delta_a).norm()<<std::endl;
    }while(k<=max_iter || (delta_new-delta_a).norm() > 1.0e-3);
    update();
    std::cout<<delta_a<<std::endl;
    std::cout<<delta_b<<std::endl;
  }
  void update() {
    //modify camera_data
    for(int cid=0; cid<camera_number; cid++) {
      Camera cam;
      cam = camera_data[cid];
      Eigen::MatrixXd temp = delta_a.block<6,1>(cid*6,0);
      cam.r1 += temp(0,0);
      cam.r2 += temp(1,0);
      cam.r3 += temp(2,0);
      
      cam.t1 += temp(3,0);
      cam.t2 += temp(4,0);
      cam.t3 += temp(5,0);
      
      camera_data[cid] = cam;
    }
    for(int pid=0; pid<point_number; pid++) {
      Point3f pt;
      pt = point_data[pid];
      pt.x += delta_b.block<3,1>(pid*3,0)(0,0);
      pt.x += delta_b.block<3,1>(pid*3,0)(1,0);
      pt.x += delta_b.block<3,1>(pid*3,0)(2,0);
      point_data[pid] = pt;
    }
    
  }
  void release() {
    U.clear();
    V.clear();
    W.clear();
    EpsilonCam.clear();
    EpsilonPt.clear();
  }
  void test() {
    
    solve_lm();
    
    //std::cout<<delta_a<<std::endl;
    //std::cout<<delta_b<<std::endl;    
    
    
  }
  

};

#endif
