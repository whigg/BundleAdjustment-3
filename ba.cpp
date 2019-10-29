#include <stdio.h>
//#include "ba.hpp"
#include "ba.hpp"
#include <Eigen/Eigen>
//#include <opencv2/opencv.hpp>

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
  // Eigen::VectorXd delta_a(3);
  // Eigen::VectorXd delta_b;
  // int n = 5;
  // delta_a = Eigen::VectorXd(n);
  // for(int i=0; i<n; i++) {
  //   delta_a(i) = double(i);
  // }
  // std::cout<<delta_a<<std::endl;
  //problem.test();
  //cv::Mat src = cv::Mat_<float>(3,1);
  //src.at<float>(0,0) = 1.5741515942940262e-02;
  //src.at<float>(1,0) = -1.2790936163850642e-02;
  //src.at<float>(2,0) = -4.4008498081980789e-03;
  //cv::Mat dst;
  //cv::Rodrigues(src, dst);
  //std::cout<<dst<<std::endl;
  int pid = 2;
  int cid = 3;
  printf("at %d, %d = %lf\n",pid, cid, problem.repj_err(pid,cid));
  Eigen::Vector3d a ;
  a<<1.0, 1.0, 1.0;
  printf("norm = %lf\n", a.norm());
  Eigen::VectorXd b ;
  b = a.transpose()*a;
  printf("norm = %lf\n",sqrt(b(0,0)));
  return 0;
}
