#include <stdio.h>
//#include "ba.hpp"
#include "public.hpp"
#include <Eigen/Eigen>

int main(int argc, char **argv) {
  BAProblem problem("./data/problem-49-7776-pre.txt");
  int pid = 7000;
  int cid = 1;
  //scanf("%d,%d", &pid, &cid);
  std::cout<<problem.W(pid, cid) <<std::endl;
  return 0;
}
