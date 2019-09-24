#ifndef PUBLIC_HPP
#define PUBLIC_HPP
#include <stdio.h>
#include <vector>

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
public:
  float perr(int const pid, int const cid ) {//calculate the projection error of point pid projecting to camera(cid) 
    //TODO reprojection i,j
    
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
    }
    for(int i=0; i<point_number; i++) {
      Point3f pt;
      fscanf(fp, "%f", &(pt.x));
      fscanf(fp, "%f", &(pt.y));
      fscanf(fp, "%f", &(pt.z));
    }
    fclose(fp);
  }
};

#endif
