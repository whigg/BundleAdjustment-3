#ifndef PUBLIC_HPP
#define PUBLIC_HPP
#include <stdio.h>
#include <vector>
#include <opencv2/opencv.hpp>

#define NOT_FOUND -1
float delta = 1.0e-5;

struct Point3f{
  float x, y, z;
};

struct Camera {
  float r1, r2, r3, t1, t2, t3, f, k1,k2;
  
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
  std::vector<Point3f> point_data;
  std::vector<Camera> camera_data;
  std::vector<Observation> observation_data;
  std::vector<int> blocks;// blocks index for accelerating find projection(i,j)
  std::vector<int> base;//index of observation_data
public:
  float pjerr(int const pid, int const cid ) {//calculate the projection error of point pid projecting to camera(cid) 
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
    return d1+d2;
  }
  float pjerr(Camera const &cam, Point3f const &pt, Observation const & obs) {
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
    return d1+d2;
  }
  cv::Mat Jac(int const pid, int const cid) {
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
    cam2.r1 += delta;
    float pj1 = pjerr(cam, pt, obs);
    float pj2 = pjerr(cam2, pt, obs);
    Jcam.at<float>(0,0) = (pj2 - pj1)/delta;
    cam2.r1 = cam.r1;

    cam2.r2 += delta;
    pj1 = pjerr(cam, pt, obs);
    pj2 = pjerr(cam2, pt, obs);
    Jcam.at<float>(0,1) = (pj2 - pj1)/delta;
    cam2.r2 = cam.r2;

    cam2.r3 += delta;
    pj1 = pjerr(cam, pt, obs);
    pj2 = pjerr(cam2, pt, obs);
    Jcam.at<float>(0,2) = (pj2 - pj1)/delta;
    cam2.r3 = cam.r3;
    //calculate Jt
    cam2.t1 += delta;
    pj1 = pjerr(cam, pt, obs);
    pj2 = pjerr(cam2, pt, obs);
    Jcam.at<float>(0,3) = (pj2 - pj1)/delta;
    cam2.t1 = cam.t1;

    cam2.t2 += delta;
    pj1 = pjerr(cam, pt, obs);
    pj2 = pjerr(cam2, pt, obs);
    Jcam.at<float>(0,4) = (pj2 - pj1)/delta;
    cam2.t2 = cam.t2;

    cam2.t3 += delta;
    pj1 = pjerr(cam, pt, obs);
    pj2 = pjerr(cam2, pt, obs);
    Jcam.at<float>(0,5) = (pj2 - pj1)/delta;
    cam2.t3 = cam.t3;
    return Jcam;
  } 
  cv::Mat Jap(int const pid, int const cid) {
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
    pt2.x += delta;
    float pj1 = pjerr(cam, pt, obs);
    float pj2 = pjerr(cam, pt2, obs);
    Jpt.at<float>(0,0) = (pj2 - pj1)/delta;
    pt2.x = pt.x;

    pt2 = pt;
    pt2.y += delta;
    pj1 = pjerr(cam, pt, obs);
    pj2 = pjerr(cam, pt2, obs);
    Jpt.at<float>(0,1) = (pj2 - pj1)/delta;
    pt2.y = pt.y;

    pt2 = pt;
    pt2.z += delta;
    pj1 = pjerr(cam, pt, obs);
    pj2 = pjerr(cam, pt2, obs);
    Jpt.at<float>(0,2) = (pj2 - pj1)/delta;
    pt2.z = pt.z;
    //std::cout<<Jpt<<std::endl;
    return Jpt;
    
  }
  void Hcam(int const pid, cv::Mat & hes) {//calculate U_i
      std::vector<int> cid;
      cid.clear();
      this->at(pid, cid);
      hes = cv::Mat();
  }
  BAProblem(char const * file_name ) {//load data
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
    for(int i=0; i<observation_number; i++) {
      Observation obs;
      fscanf(fp, "%d", &(obs.cid));
      fscanf(fp, "%d", &(obs.pid));
      fscanf(fp, "%f", &(obs.x));
      fscanf(fp, "%f", &(obs.y));
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
    for(int i=0; i<blocks.size(); i++) {
      base.push_back(base[i] + blocks[i]);
    }
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
      camera_data.push_back(cam);
    } 
    for(int i=0; i<point_number; i++) {
      Point3f pt;
      fscanf(fp, "%f", &(pt.x));
      fscanf(fp, "%f", &(pt.y));
      fscanf(fp, "%f", &(pt.z));
      point_data.push_back(pt);
    }
    fclose(fp);
  }
  int at(int const pid, int const cid) {
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
  void at(int const pid, std::vector<int> &cid) {
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
  void print_blocks() {
    //printf("%d\n", base.size());
    // for(int i=0; i<blocks.size(); i++) {
    //   printf("%d %d, ", base[i],blocks[i]);
    // }
    // printf("\n");
  }
};
  


#endif
