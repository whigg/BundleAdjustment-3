#include <stdio.h>
//#include "ba.hpp"
#include "public.hpp"
#include <Eigen/Eigen>

int main(int argc, char **argv) {
  BAProblem problem("./data/problem-49-7776-pre.txt");
  //int pid = 7000;
  int cid = 1;
  cv::Mat epsa = problem.epsilon_cam(cid);
  std::cout<<epsa<<std::endl;
  return 0;
}
