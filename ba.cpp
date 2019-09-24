#include <stdio.h>
#include "ba.hpp"
#include "public.hpp"
#include <Eigen/Eigen>


int main(int argc, char **argv) {
    //Function func_arr[3];
    std::vector<unsigned> cameraIndex;
    std::vector<unsigned> pointIndex;
    cv::Mat observation;
    cv::Mat camera;
    cv::Mat points;
    LoadData("./data/problem-49-7776-pre.txt",
	     cameraIndex, pointIndex, observation, camera, points);
    
    printf("size point index = %lu\n", pointIndex.size());
    printf("size camera index = %lu\n", cameraIndex.size());
    printf("size of camera = (%u, %u)\n", camera.rows, camera.cols);
    printf("size of points = (%u, %u)\n", points.rows, points.cols);
    printf("size of observation = (%u, %u)\n", observation.rows, observation.cols);
    double total_cost = TotalCost(reprojection_err, camera, points, cameraIndex, pointIndex,
				  observation);
    printf("total cost = %lf\n",total_cost);
    cv::Mat JacMatrix = Jacobian(reprojection_err, camera, points, cameraIndex, pointIndex,
				 observation);
    printf("Jacobian column size = (%d, %d)\n", JacMatrix.rows, JacMatrix.cols);
    
    cv::namedWindow("bitmap", cv::WINDOW_NORMAL);
    cv::imshow("bitmap", JacMatrix);
    //cv::imwrite("./bitmap.jpg", JacMatrix);
    cv::waitKey(0);
    return 0;
}
