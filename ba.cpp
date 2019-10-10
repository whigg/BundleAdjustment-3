#include <stdio.h>
//#include "ba.hpp"
#include "public.hpp"
#include <Eigen/Eigen>

int main(int argc, char **argv) {
  BAProblem problem("./data/problem-49-7776-pre.txt");
  
  // cv::Mat sm = problem.SM();
  // cv::Mat em = problem.em();
  // std::cout<<sm.rows<<", "<<sm.cols<<std::endl;
  // cv::Mat delta_a;
  // cv::solve(sm, em, delta_a );

  //std::cout<<delta_a<<std::endl;
  problem.solve_lineq();
  return 0;
}
