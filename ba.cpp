#include <stdio.h>
#include "ba.hpp"
#include "public.hpp"
#include <Eigen/Eigen>


int main(int argc, char **argv) {
  BAProblem problem("./data/problem-49-7776-pre.txt");
  problem.print_blocks();
  return 0;
}
