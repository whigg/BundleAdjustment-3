#ifndef BA_HPP
#define BA_HPP

#include <vector>
#include <stdio.h>
#include <opencv2/opencv.hpp>


#define INFINITESIMAL 1.0e-6
#define N_CAMERA_PARAMETER 9
#define N_POINT_3D         3
#define N_OBSERVED_2D      2


typedef double (*Function) (std::vector<double> &, std::vector<double> &,
                            std::vector<double> &);

double reprojection_err(std::vector<double> &camera, std::vector<double> &point_3d, std::vector<double> & observed_2d) {
    //TODO ..
    
    if(camera.size() != 9) {
        printf("error camera size\n");
        exit(-1);
    }
    cv::Mat rvec = cv::Mat_<double>(1,3);
    rvec.at<double>(0,0) = camera[0];
    rvec.at<double>(0,1) = camera[1];
    rvec.at<double>(0,2) = camera[2];
    cv::Mat t = cv::Mat_<double>(3,1);
    t.at<double>(0,0) = camera[3];
    t.at<double>(1,0) = camera[4];
    t.at<double>(2,0) = camera[5];
    double f = camera[6];
    double k1 = camera[7];
    double k2 = camera[8];
    cv::Mat rotm = cv::Mat_<double>(3,3);
    cv::Rodrigues(rvec, rotm);//rodrigues transformation
    cv::Mat Xw = cv::Mat_<double>(3,1);
    Xw.at<double>(0,0) = point_3d[0];
    Xw.at<double>(1,0) = point_3d[1];
    Xw.at<double>(2,0) = point_3d[2];
    cv::Mat Xc = rotm * Xw + t;
    Xc = -Xc/Xc.at<double>(2,0);

    double l1 = cv::norm(Xc, cv::NORM_L2);
    double l2 = l1*l1;
    double rho = 1.0 + k1 * l1 + k2 * l2;
    cv::Mat Xp = f * rho * Xc;
    double d1, d2;
    d1 = Xp.at<double>(0,0) - observed_2d[0];
    d2 = Xp.at<double>(1,0) - observed_2d[1];
    return sqrt(d1*d1 + d2*d2);
}
/*******************************************************************************************************************************/
/* camera          : camera parameters to be optimized n*9 dimension matrix                                                     *
*  point_3d        : 3D points parameters to be optimized, n*3 dimension matrix                                                 *
*  observed_2d_data: Observed 2D data with dimension of n*2                                                                     *
*  cameraIndex     : Camera index, from 0 ... to number of cameras -1, which == camera.rows()-1                                 *
*  return value    : n * (9*nCam + 3 * nPoints) Jacobian matrix where nCam == number of camera, nPoints == number of 3D points  *
*********************************************************************************************************************************/
cv::Mat Jacobian(Function, cv::Mat &cameras, cv::Mat & point_3ds,
                            std::vector<unsigned> & cameraIndex, std::vector<unsigned> & pointIndex, 
                            cv::Mat & observed_2d) {
    //TODO ...
    cv::Mat retval;
    return retval;
}
double TotalCost(Function, cv::Mat &cameras, cv::Mat & point_3ds,
                            std::vector<unsigned> & cameraIndex, std::vector<unsigned> & pointIndex, 
                            cv::Mat & observed_2d) {
    unsigned nItem = cameraIndex.size();
    
    double sum = 0.0;
    std::vector<double> cam;
    std::vector<double> pt_3d;
    std::vector<double> obs_2d;
    for(unsigned i = 0; i< nItem; i++ ) {
        cam.clear();
        pt_3d.clear();
        obs_2d.clear();
        //get camera parameters for camera[cam_idx]
        for(int j=0; j<N_CAMERA_PARAMETER; j++) {
            cam.push_back(cameras.at<double>(cameraIndex[i], j));
        }
        //get point 3d parameters at index i
        for(int j=0; j<N_POINT_3D; j++) {
            pt_3d.push_back(point_3ds.at<double>(pointIndex[i], j));
        }
        for(int j=0; j<N_OBSERVED_2D; j++) {
            obs_2d.push_back(observed_2d.at<double>(i,j));
        }
        double single_error = reprojection_err(cam, pt_3d, obs_2d);
        //printf("single error =%lf\n", single_error);
        sum += single_error;
    }
}

void LoadData(char const *file, std::vector<unsigned> &camIndex, std::vector<unsigned> & pointIndex, 
                cv::Mat & observe, cv::Mat &camera, cv::Mat &points) {
    FILE *fp = fopen(file, "r");
    if(!fp) {
        perror("error reading file\n");
        exit(-1);
    }
    int nCamera;//number of cameras
    int nPoint; //number of 3d points
    int nObserve; //number of observations
    /*read first line*/
    fscanf(fp, "%d", &nCamera);
    fscanf(fp, "%d", &nPoint);
    fscanf(fp, "%d", &nObserve);
    printf("nCamera = %d\n", nCamera);
    observe = cv::Mat_<double>(nObserve, 2);//n*2 matrix, n is number of observations
    camIndex.clear();
    pointIndex.clear();
    //read the second part, camera index, point index, observation x and y
    for(int i=0;i<nObserve; i++) {
        unsigned cidx;
        unsigned pidx;
        fscanf(fp, "%u",&cidx );
        fscanf(fp, "%u",  &pidx);
        camIndex.push_back(cidx);
        pointIndex.push_back(pidx);
        double x,y;
        fscanf(fp, "%lf",&x);
        fscanf(fp, "%lf",&y);
        observe.at<double>(i, 0) = x;
        observe.at<double>(i, 1) = y;
    }
    //read the third part , camera parameters
    camera = cv::Mat_<double>(nCamera, 9);//n*9 dimension
    for(int i=0; i<nCamera; i++) {
        for(int j=0; j<9; j++) {
            double temp = 0.0;
            fscanf(fp, "%lf", &temp);
            camera.at<double>(i,j) = temp;
        }
    }
    //read the forth part, point x,y,z
    points = cv::Mat_<double>(nPoint, 3);//n*3
    for(int i=0; i<nPoint; i++ ) {
        for(int j=0; j<3; j++) {
            double temp = 0.0;
            fscanf(fp, "%lf", &temp);
            points.at<double>(i,j) = temp;
        }
    }
    fclose(fp);
}


#endif
