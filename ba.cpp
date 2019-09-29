#include <stdio.h>
//#include "ba.hpp"
#include "public.hpp"
#include <Eigen/Eigen>

int main(int argc, char **argv) {
  BAProblem problem("./data/problem-49-7776-pre.txt");
  cv::Mat U = problem.U(1);
  cv::Mat Us = problem.Us(1);
  float factor = 1e-6;
  Us = factor*Us;
  std::cout<<"Us = \n"<<Us<<std::endl;
  std::cout<<"Us.inv.inv = \n"<<Us.inv(cv::DECOMP_SVD).inv()-Us<<std::endl;
  
  return 0;
}
