#include <stdio.h>
//#include "ba.hpp"
#include "public.hpp"
#include <Eigen/Eigen>

class Cell{
  
};

int main(int argc, char **argv) {
  BAProblem problem("./data/problem-49-7776-pre.txt");
  
  // cv::Mat sm = problem.SM();
  // cv::Mat em = problem.em();
  // std::cout<<sm.rows<<", "<<sm.cols<<std::endl;
  // cv::Mat delta_a;
  //cv::solve(sm, em, delta_a );

  //std::cout<<delta_a<<std::endl;
  //problem.solve_lineq();
  //problem.test();
  Eigen::VectorXd delta_a(3);
  Eigen::VectorXd delta_b;
  int n = 5;
  delta_a = Eigen::VectorXd(n);
  for(int i=0; i<n; i++) {
    delta_a(i) = double(i);
  }
  std::cout<<delta_a<<std::endl;
  return 0;
}
