#include <stdio.h>
//#include "ba.hpp"
#include "public.hpp"
#include <Eigen/Eigen>

int main(int argc, char **argv) {
  BAProblem problem("./data/problem-49-7776-pre.txt");
  
  cv::Mat sm = problem.SM();
  
  std::cout<<sm.size()<<std::endl;
  //std::cout<<problem.S(1,1)<<std::endl;
  return 0;
}
