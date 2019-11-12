#include <stdio.h>
#include "ba.hpp"
#include <Eigen/Eigen>

class Cell{
  
};

int main(int argc, char **argv) {
  BAProblem problem("./data/problem-49-7776-pre.txt");
  
  problem.test();
  //Eigen::MatrixXd a;
  Eigen::Matrix<double,3,4> a;
  a << 1, 2, 3, 4, 4, 3, 2, 1, 4,5, 2,1;
  printf("maximum = %lf\n", problem.maxd(a));

  return 0;
}
