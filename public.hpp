#ifndef PUBLIC_HPP
#define PUBLIC_HPP
#include <stdio.h>

#include <opencv2/opencv.hpp>
cv::Mat selfMultiply(cv::Mat const & m ) {/// return a' * a
    int row = m.rows;
    int col = m.cols;
    /*
    分块矩阵[a,b
            c,d]
    */
    //cv::Mat a,b,c,d;
    int rowa,rowd,cola,cold;
    rowa = row/2;
    cola = col/2;
    rowd = row - rowa;
    cold = col - cold;
    cv::Mat a = cv::Mat_<double>(rowa, cola);
    cv::Mat b = cv::Mat_<double>(rowa, cold);
    cv::Mat c = cv::Mat_<double>(rowd, cola);
    cv::Mat d = cv::Mat_<double>(rowd, cold);
    //copy m to [a,b;c,d]

}

#endif