#ifndef PUBLIC_HPP
#define PUBLIC_HPP
#include <stdio.h>
#include <vector>
#include <iterator>
#include <opencv2/opencv.hpp>

#define NOT_FOUND -1

typedef std::vector<std::vector<int> > Table;
double delta = 0.00002;
//typedef float CameraPose[6]; //camera Pose parameters r1, r2, r3, t1, t2, t3, to be optimized;
//typedef float PointParameter[3];//3D point x, y, z to be optimized

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
struct BAProblem {
private:
  int camera_number;
  int point_number;
  int observation_number;
  float tau; // if x* is a good approximation tau can be set as a small float number like 1e-6, otherwise
               // set tau = 1e-3 or even 1
  std::vector<Point3f> point_data;
  std::vector<Camera> camera_data;
  std::vector<Observation> observation_data;
  std::vector<cv::Mat> da;//delta a 
  std::vector<cv::Mat> db;//delta b
  cv::Mat p;
  cv::Mat g;
  cv::Mat dp;
  int k;
  int v ;
  std::vector<int> blocks;// blocks index for accelerating find projection(i,j)
  std::vector<int> base;//index of observation_data
  Table captt; // camera point table captt[i] return vector containing points projected to camera i
  //cv::Mat pcam;
  //cv::Mat ppt; 
public:
  float pjerr(int const pid, int const cid ) const {//calculate the projection error of point pid projecting to camera(cid) 
    //TODO reprojection i,j
    int obs_i = this->at(pid, cid);
    if(obs_i == NOT_FOUND) {
      return 0.0;
    }
    float obsx, obsy;
    obsx = observation_data[obs_i].x;
    obsy = observation_data[obs_i].y;
    
    Camera cam = camera_data[cid];
    Point3f pt = point_data[pid];
    //calculate projection error
    cv::Mat rvec = cv::Mat_<float>(1,3);
    rvec.at<float>(0,0) = cam.r1;
    rvec.at<float>(0,1) = cam.r2;
    rvec.at<float>(0,2) = cam.r3;
    cv::Mat t = cv::Mat_<float>(3,1);
    t.at<float>(0,0) = cam.t1;
    t.at<float>(0,1) = cam.t2;
    t.at<float>(0,2) = cam.t3;
    float f = cam.f;
    float k1 = cam.k1;
    float k2 = cam.k2;
    cv::Mat rotm = cv::Mat_<float>(3,3);
    cv::Rodrigues(rvec, rotm);//rodrigues transformation
    cv::Mat Xw = cv::Mat_<float>(3,1);
    Xw.at<float>(0,0) = pt.x;
    Xw.at<float>(1,0) = pt.y;
    Xw.at<float>(2,0) = pt.z;
    cv::Mat Xc = rotm * Xw + t;
    Xc = -Xc / Xc.at<float>(2,0);
    float l1 = cv::norm(Xc, cv::NORM_L2);
    float l2 = l1*l1;
    float rho = 1.0 + k1 * l1*l1 + k2 * l2*l2;
    cv::Mat Xp = f * rho * Xc;
    float d1 = fabs(Xp.at<float>(0,0) - obsx);
    float d2 = fabs(Xp.at<float>(1,0) - obsy);
    return sqrt(d1*d1+d2*d2);
  }
  float pjerr(Camera const &cam, Point3f const &pt, Observation const & obs) const {
    cv::Mat rvec = cv::Mat_<float>(1,3);
    rvec.at<float>(0,0) = cam.r1;
    rvec.at<float>(0,1) = cam.r2;
    rvec.at<float>(0,2) = cam.r3;
    cv::Mat t = cv::Mat_<float>(3,1);
    t.at<float>(0,0) = cam.t1;
    t.at<float>(0,1) = cam.t2;
    t.at<float>(0,2) = cam.t3;
    float f = cam.f;
    float k1 = cam.k1;
    float k2 = cam.k2;
    cv::Mat rotm = cv::Mat_<float>(3,3);
    cv::Rodrigues(rvec, rotm);//rodrigues transformation
    cv::Mat Xw = cv::Mat_<float>(3,1);
    Xw.at<float>(0,0) = pt.x;
    Xw.at<float>(1,0) = pt.y;
    Xw.at<float>(2,0) = pt.z;
    cv::Mat Xc = rotm * Xw + t;
    Xc = -Xc / Xc.at<float>(2,0);
    float l1 = cv::norm(Xc, cv::NORM_L2);
    float l2 = l1*l1;
    float rho = 1.0 + k1 * l1*l1 + k2 * l2*l2;
    cv::Mat Xp = f * rho * Xc;
    float d1 = fabs(Xp.at<float>(0,0) - obs.x);
    float d2 = fabs(Xp.at<float>(1,0) - obs.y);
    return sqrt(d1*d1+d2*d2);
  }

  cv::Mat Jac(int const pid, int const cid) const {
    if(pid >= point_number || cid >= camera_number) {
      printf("index exceeded\n");
      exit(-1);
    }
    cv::Mat Jcam;
    Camera cam;
    Point3f pt;
    Observation obs;
    cam = camera_data[cid];
    pt = point_data[pid];
    int obs_i = this->at(pid, cid);
    Jcam = cv::Mat::zeros(1,6,CV_32F);
    if(obs_i == NOT_FOUND) {
      return Jcam;
    }
    obs.pid = pid;
    obs.cid = cid;
    obs.x = observation_data[obs_i].x;
    obs.y = observation_data[obs_i].y;
    //jcam 
    //calculate Jr
    Camera cam2 = cam;
    Camera cam1 = cam;
    cam2.r1 += delta;
    cam1.r1 -= delta;
    float pj1 = pjerr(cam1, pt, obs);
    float pj2 = pjerr(cam2, pt, obs);
    Jcam.at<float>(0,0) = (pj2 - pj1)/(delta*2);
    cam2.r1 = cam.r1;
    cam1.r1 = cam.r1;

    cam2.r2 += delta;
    cam1.r2 -= delta;
    pj1 = pjerr(cam1, pt, obs);
    pj2 = pjerr(cam2, pt, obs);
    Jcam.at<float>(0,1) = (pj2 - pj1)/(2*delta);
    cam2.r2 = cam.r2;
    cam1.r2 = cam.r2;

    cam2.r3 += delta;
    cam1.r3 -= delta;

    pj1 = pjerr(cam1, pt, obs);
    pj2 = pjerr(cam2, pt, obs);
    Jcam.at<float>(0,2) = (pj2 - pj1)/(2*delta);
    cam2.r3 = cam.r3;
    cam1.r3 = cam.r3;

    //calculate Jt
    cam2.t1 += delta;
    cam1.t1 -= delta;

    pj1 = pjerr(cam1, pt, obs);
    pj2 = pjerr(cam2, pt, obs);
    Jcam.at<float>(0,3) = (pj2 - pj1)/(2*delta);
    cam2.t1 = cam.t1;
    cam1.t1 = cam.t1;

    cam2.t2 += delta;
    cam1.t2 -= delta;

    pj1 = pjerr(cam1, pt, obs);
    pj2 = pjerr(cam2, pt, obs);
    Jcam.at<float>(0,4) = (pj2 - pj1)/(2*delta);
    cam2.t2 = cam.t2;
    cam1.t2 = cam.t2;

    cam2.t3 += delta;
    cam1.t3 -= delta;
    pj1 = pjerr(cam1, pt, obs);
    pj2 = pjerr(cam2, pt, obs);
    Jcam.at<float>(0,5) = (pj2 - pj1)/(2*delta);
    cam2.t3 = cam.t3;
    cam1.t3 = cam.t3;
    return Jcam;
  } 
  cv::Mat Jap(int const pid, int const cid) const {
    if(pid >= point_number || cid >= camera_number) {
      printf("index exceeded\n");
      exit(-1);
    }
    cv::Mat Jpt;
    Camera cam;
    Point3f pt;
    Observation obs;
    cam = camera_data[cid];
    pt = point_data[pid];
    int obs_i = this->at(pid, cid);
    Jpt = cv::Mat::zeros(1,3, CV_32F);
    if(obs_i == NOT_FOUND) {
      return Jpt;
    }
    obs.pid = pid;
    obs.cid = cid;
    obs.x = observation_data[obs_i].x;
    obs.y = observation_data[obs_i].y;
    //jpt
    Point3f pt2 = pt;
    Point3f pt1 = pt;
    pt2.x += delta;
    pt1.x -= delta;
    float pj1 = pjerr(cam, pt1, obs);
    float pj2 = pjerr(cam, pt2, obs);
    Jpt.at<float>(0,0) = (pj2 - pj1)/(2*delta);
    pt2.x = pt.x;
    pt1.x = pt.x;

    pt2 = pt;
    pt1 = pt;
    pt2.y += delta;
    pt1.y -= delta;
    pj1 = pjerr(cam, pt1, obs);
    pj2 = pjerr(cam, pt2, obs);
    Jpt.at<float>(0,1) = (pj2 - pj1)/(2*delta);
    pt2.y = pt.y;
    pt1.y = pt.y;

    pt2 = pt;
    pt1 = pt;
    pt2.z += delta;
    pt1.z -= delta;

    pj1 = pjerr(cam, pt, obs);
    pj2 = pjerr(cam, pt2, obs);
    Jpt.at<float>(0,2) = (pj2 - pj1)/(2*delta);
    pt2.z = pt.z;
    pt1.z = pt.z;
    //std::cout<<Jpt<<std::endl;
    return Jpt;
    
  }
  cv::Mat U(int const cid) const {//calculate U_j
      if(cid >= camera_number) {printf("index exceeded\n"); exit(-1);}
      std::vector<int> pid_vec;
      pid_vec.clear();
      pid_vec = captt[cid];
      cv::Mat Uj = cv::Mat::zeros(6,6,CV_32F);
      for(std::vector<int>::iterator it = pid_vec.begin(); it!= pid_vec.end(); it++) {
        int pid = *it;
        cv::Mat jac = this->Jac(pid, cid);
        Uj += jac.t() * jac;
      }
      return Uj;
  }
  cv::Mat Us(int const cid) const { // caluclate U_j^*
    // TODO
    cv::Mat Uj = U(cid);
    double u = tau * maxd(Uj);
    //printf("u = %f\n", u);
    cv::Mat mu = u * cv::Mat::eye(Uj.rows,Uj.cols, CV_32F);
    return mu+Uj;
  }
  float maxd(cv::Mat const diag) const {
    int m = diag.rows <= diag.cols ?diag.rows : diag.cols;
    double maximum = diag.at<float>(0,0);
    for(int i=0; i<m; i++){
      maximum = maximum > diag.at<float>(i,i) ? maximum:diag.at<float>(i,i);
    }
    return maximum;
  }
  cv::Mat V(int const pid) const {//calculate V_i
    if(pid >= point_number) {printf("index exceeded\n"); exit(-1);}
    std::vector<int> vcid;
    vcid.clear();
    this->at(pid, vcid);
    //printf("vcid size = %lu\n", vcid.size());
    cv::Mat Vi = cv::Mat::zeros(3,3,CV_32F) ;
    for(std::vector<int>::iterator iter = vcid.begin(); iter != vcid.end(); iter++) {
      int cid = *iter;
      //printf("cid = %d\n", cid);
      cv::Mat Bij = this->Jap(pid, cid);
      Vi += Bij.t() * Bij;
    }
    return Vi;
  }
  cv::Mat Vs(int const pid) const {//calculate Vi*
    cv::Mat Vi = V(pid);
    double u = tau * maxd(Vi);
    cv::Mat mu = u * cv::Mat::eye(Vi.rows, Vi.cols, CV_32F);
    return mu + Vi;
  }
  cv::Mat Y(int const pid, int const cid) const {
    cv::Mat Vsi = this->Vs(pid);
    cv::Mat wij = this->W(pid, cid);
    return wij * Vsi.inv();
  }
  cv::Mat W(int const pid, int const cid) const {
    if(pid >= point_number || cid >= camera_number) {
      printf("index exceeded\n");
      exit(-1);
    }
    cv::Mat wij = cv::Mat::zeros(6,3, CV_32F);
    if(this->at(pid,cid)== NOT_FOUND) 
      return wij;
    cv::Mat Aij = Jac(pid, cid);
    cv::Mat Bij = Jap(pid, cid);
    wij = Aij.t() * Bij;
    return wij;
  }
  cv::Mat epsilon_cam(int const cid) const {
    //find pid correspoinding to cid
    std::vector<int> pid_vec;
    pid_vec.clear();
    pid_vec = this->captt[cid];
    cv::Mat epsa = cv::Mat::zeros(6,1,CV_32F);
    for(auto it = pid_vec.begin(); it != pid_vec.end(); it++) {
      int pid = *it;
      cv::Mat jacm = this->Jac(pid, cid);
      float e = this->pjerr(pid, cid);
      epsa += jacm.t() * e;
    }
    return epsa;
  }
  cv::Mat epsilon_pt(int const pid) const {
    std::vector<int> cid_vec;
    cid_vec.clear();
    this->at(pid, cid_vec);
    cv::Mat epsb = cv::Mat::zeros(3,1,CV_32F);
    for(auto it = cid_vec.begin(); it!=cid_vec.end(); it++) {
      int cid = *it;
      cv::Mat japm = this->Jap(pid, cid);
      float e = this->pjerr(pid, cid);
      epsb += japm.t() * e;
    }
    return epsb;
  }
  cv::Mat S(int const j, int const k) const {
    //TODO
    cv::Mat Usjk = Us(j);//j==k
    int r = Usjk.rows;
    int c = Usjk.cols;
    if(j != k) {
      Usjk = cv::Mat::zeros(r,c,CV_32F);
    }
    cv::Mat YW = cv::Mat::zeros(r, c, CV_32F);
    int cid = j;
    std::vector<int> vpid = captt[cid];
    for(std::vector<int>::iterator it=vpid.begin(); it!=vpid.end(); it++) {
      int i = *it;
      YW += Y(i,j) * W(i,k).t();
    }
    return Usjk - YW;
  }
  cv::Mat e(int const cid)const {
    //epsilon_cam(cid) - sum(Y(pid, cid) * epsilon_pt(pid));
    cv::Mat s= cv::Mat::zeros(6,1, CV_32F) ;
    std::vector<int> vpid = captt[cid];
    for(std::vector<int>::iterator it = vpid.begin(); it!= vpid.end(); it++) {
      int pid = *it;
      s += Y(pid, cid)*epsilon_pt(pid);
    }
    return epsilon_cam(cid) - s;
  }
  cv::Mat em() const {
    std::vector<cv::Mat> ve;
    for(int cid=0; cid<camera_number; cid++) {
      ve.push_back(e(cid));
    }
    cv::Mat re;//return e
    cv::vconcat(ve, re);
    return re;
  }
  cv::Mat SM() {
    cv::Mat sm;
    std::vector<cv::Mat> col;
    std::vector<cv::Mat> row;
    row.clear();
    for(int i=0; i<camera_number; i++ ) {
      col.clear();
      cv::Mat colm;
      for(int j=0; j<camera_number; j++) {
        cv::Mat s = S(i,j);
        col.push_back(s);
      }
      cv::hconcat(col, colm);
      row.push_back(colm);
    }
    cv::vconcat(row, sm);
    return sm;
  }
  void solve_lineq(){//solve linear equation and store delta a and delta b to da,db
    cv::Mat sm = SM();
    //cv::Mat sm = cv::Mat_<float>(294,294);
    cv::Mat em = this->em();
    cv::Mat a;
    cv::solve(sm, em, a);
    std::cout<<"a = \n"<<a<<std::endl;
    std::vector<cv::Mat> dta;//delta a
    dta.clear();
    cv::Mat temp = cv::Mat_<float>(6,1);
    for(int i=0; i<camera_number; i++) {
      for(int j=0; j<6; j++) {
        temp.at<float>(j,0) = a.at<float>(i*camera_number+j,0);
      }
      dta.push_back(temp.clone());
    }
    //std::cout<<"debug"<<std::endl;
    cv::Mat b;
    std::vector<cv::Mat> dtb;
    temp = cv::Mat_<float>(3,1);
    for(int pid = 0; pid < point_number; pid++) {
      std::vector<int> vcid;
      this->at(pid, vcid);
      cv::Mat vsinv = this->Vs(pid).inv();
      cv::Mat epsb = this->epsilon_pt(pid);
      cv::Mat sum = cv::Mat::zeros(3,1,CV_32F);
      //std::cout<<W(pid, 0).rows<<", "<< W(pid, 0).cols<<std::endl;
      for(auto it = vcid.begin(); it != vcid.end(); it++) {
        int cid = *it;

        sum += this->W(pid, cid).t()*dta[cid];
      }
      std::cout<<"sum = \n";
      std::cout<<sum<<std::endl;
      dtb.push_back(vsinv*(epsb - sum));
    }
    
    cv::vconcat(dtb, b);
    cv::vconcat(a,b, dp);
  }
  void solve_lm(){ //solve linear
    //TODO
    return;
  }
  BAProblem(char const * file_name ) {//load data
    tau = 1.0e-6; // set tau here
    FILE * fp = fopen(file_name, "r");
    if(!fp) {
      perror("error opening file\n");
      exit(-1);
    }
    fscanf(fp, "%d", &camera_number);
    fscanf(fp, "%d", &point_number);
    fscanf(fp, "%d", &observation_number);
    camera_data.clear();
    point_data.clear();
    observation_data.clear();
    blocks.clear();
    base.clear();
    int current_pid = 0;
    int current_pid_num = 0;
    for(int i=0; i<camera_number; i++) {
      std::vector<int> empty_vec;
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
	      //printf("obs.pid = %d\n", obs.pid);
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
      base.push_back(base[i] + blocks[i]);
    }
    this->da.clear();
    for(int i=0; i<camera_number; i++) {
      Camera cam;
      fscanf(fp, "%f", &(cam.r1));
      fscanf(fp, "%f", &(cam.r2));
      fscanf(fp, "%f", &(cam.r3));
      fscanf(fp, "%f", &(cam.t1));
      fscanf(fp, "%f", &(cam.t2));
      fscanf(fp, "%f", &(cam.t3));
      fscanf(fp, "%f", &(cam.f));
      fscanf(fp, "%f", &(cam.k1));
      fscanf(fp, "%f", &(cam.k2));
      cv::Mat cp = cv::Mat_<float>(6,1);
      cp.at<float>(0,0) = cam.r1;
      cp.at<float>(1,0) = cam.r2;
      cp.at<float>(2,0) = cam.r3;
      cp.at<float>(3,0) = cam.t1;
      cp.at<float>(4,0) = cam.t2;
      cp.at<float>(5,0) = cam.t3;
      da.push_back(cp);
      camera_data.push_back(cam);
    } 
    this->db.clear();
    for(int i=0; i<point_number; i++) {
      Point3f pt;
      cv::Mat pp = cv::Mat_<float>(3,1);
      fscanf(fp, "%f", &(pt.x));
      fscanf(fp, "%f", &(pt.y));
      fscanf(fp, "%f", &(pt.z));
      point_data.push_back(pt);
      pp.at<float>(0,0) = pt.x;
      pp.at<float>(1,0) = pt.y;
      pp.at<float>(2,0) = pt.z;
      db.push_back(pp);
    }
    fclose(fp);
  }
  int at(int const pid, int const cid) const {
    if(pid >= point_number || cid >= camera_number) {
      printf("point index or camera index exceeded\n");
      exit(-1);
    }
    int base = this->base[pid];
    int limit = this->base[pid+1];
    for(int i=base;i < limit; i++) {
      Observation obs;
      obs = observation_data[i];
      if(pid == obs.pid && cid == obs.cid) {
        return i;
      }
    }
    return NOT_FOUND;
  }
  void at(int const pid, std::vector<int> &cid) const {
    if(pid >= point_number) {
      printf("point index exceeded\n");
      exit(-1);
    }
    int base = this->base[pid];
    int limit = this->base[pid+1];
    cid.clear();
    for(int i = base; i<limit; i++) {
      Observation obs;
      obs = observation_data[i];
      cid.push_back(obs.cid);
    }
  }
  void test() {
    // cv::Mat sm = this->SM();
    // //cv::Mat sm = cv::Mat_<float>(294,294);
    // cv::Mat em = this->em();
    // cv::Mat a;
    // cv::solve(sm, em, a);
    // std::vector<cv::Mat> dta;//delta a
    // dta.clear();
    // cv::Mat temp = cv::Mat_<float>(6,1);
    // for(int i=0; i<camera_number; i++) {
    //   for(int j=0; j<6; j++) {
    //     temp.at<float>(j,0) = a.at<float>(i*camera_number+j,0);
    //   }
    //   dta.push_back(temp.clone());
    // }
    // //std::cout<<"debug"<<std::endl;
    // cv::Mat b;
    // std::vector<cv::Mat> dtb;
    // temp = cv::Mat_<float>(3,1);
    // for(int pid = 0; pid < point_number; pid++) {
    //   std::vector<int> vcid;
    //   this->at(pid, vcid);
    //   cv::Mat vsinv = this->Vs(pid).inv();
    //   cv::Mat epsb = this->epsilon_pt(pid);
    //   cv::Mat sum = cv::Mat::zeros(3,1,CV_32F);
    //   //std::cout<<W(pid, 0).rows<<", "<< W(pid, 0).cols<<std::endl;
    //   for(auto it = vcid.begin(); it != vcid.end(); it++) {
    //     int cid = *it;

    //     sum += this->W(pid, cid).t()*dta[cid];
    //   }
    //   std::cout<<"sum = \n";
    //   std::cout<<sum<<std::endl;
    //   dtb.push_back(vsinv*(epsb - sum));
    // }
    
    // cv::vconcat(dtb, b);
    // cv::vconcat(a,b, dp);
    std::cout<<this->W(0,0).t()<<std::endl;
  }
    
};
  


#endif
