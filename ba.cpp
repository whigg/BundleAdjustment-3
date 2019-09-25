#include <stdio.h>
#include "ba.hpp"
#include "public.hpp"
#include <Eigen/Eigen>


int main(int argc, char **argv) {
  BAProblem problem("./data/problem-49-7776-pre.txt");
  int pid = 0;
  int cid = 0;
  //scanf("%d,%d", &pid, &cid);
  printf("item %d,%d = %d\n",pid, cid, problem.at(pid,cid));
  printf("projection error = %f\n", problem.pjerr(pid, cid));
  printf("Jcam = \n");
  cv::Mat Jcam;
  Jcam = problem.Jac(pid, cid);
  printf("%f %f %f %f %f %f\n", Jcam.at<float>(0,0),Jcam.at<float>(0,1),
          Jcam.at<float>(0,2),Jcam.at<float>(0,3),Jcam.at<float>(0,4),
          Jcam.at<float>(0,5) );
  cv::Mat Jpt;
  printf("Jap = \n");
  Jpt = problem.Jap(pid, cid);
  printf("%f %f %f\n", Jpt.at<float>(0,0), Jpt.at<float>(0,1), Jpt.at<float>(0,2));
  return 0;
}
