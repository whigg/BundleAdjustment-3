#ifndef PUBLIC_HPP
#define PUBLIC_HPP
#include <stdio.h>

#include <opencv2/opencv.hpp>
cv::Mat selfMultiply(cv::Mat const & m ) {/// return a' * a
    //assert(mul.cols == m.rows);
    cv::Mat result = cv::Mat::zeros(m.cols, m.cols, CV_32F );
    printf("...\n");
    for(int i=0; i<m.cols; i++) {
        for(int j = 0; j<m.cols; j++) {
            for(int k = 0; k<m.rows; k++) {
                result.at<float>(i,j) += m.at<float>(k,i) * m.at<float>(k,j);//m.at<double>(k,i)
            }
        }
    }
    return result;
}
#endif