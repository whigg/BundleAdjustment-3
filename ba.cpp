#include <stdio.h>
#include "ba.hpp"


int main(int argc, char **argv) {
    //Function func_arr[3];
    std::vector<unsigned> cameraIndex;
    std::vector<unsigned> pointIndex;
    cv::Mat observation;
    cv::Mat camera;
    cv::Mat points;
    LoadData("./data/problem-49-7776-pre.txt", cameraIndex, pointIndex, observation, camera, points);
    printf("size point index = %lu\n", pointIndex.size());
    printf("size camera index = %lu\n", cameraIndex.size());
    printf("size of camera = (%u, %u)\n", camera.rows, camera.cols);
    printf("size of points = (%u, %u)\n", points.rows, points.cols);
    printf("size of observation = (%u, %u)\n", observation.rows, observation.cols);
    double total_cost = TotalCost(reprojection_err, camera, points, cameraIndex, pointIndex, observation);
    Jacobian(reprojection_err, camera, points, cameraIndex, pointIndex, observation);
    printf("total cost = %lf\n",total_cost);
    return 0;
}
