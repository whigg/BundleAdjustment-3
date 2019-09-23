#ifndef BA_HPP
#define BA_HPP

#include <vector>
#include <stdio.h>
#include <opencv2/opencv.hpp>


#define INFINITESIMAL 1.0e-6
#define N_CAMERA_PARAMETER 9
#define N_POINT_3D         3
#define N_OBSERVED_2D      2


typedef float (*Function) (std::vector<float> &, std::vector<float> &,
                            std::vector<float> &);

float reprojection_err(std::vector<float> &camera, std::vector<float> &point_3d, std::vector<float> & observed_2d) {
    
    if(camera.size() != 9) {
        printf("error camera size\n");
        exit(-1);
    }
    cv::Mat rvec = cv::Mat_<float>(1,3);
    rvec.at<float>(0,0) = camera[0];
    rvec.at<float>(0,1) = camera[1];
    rvec.at<float>(0,2) = camera[2];
    cv::Mat t = cv::Mat_<float>(3,1);
    t.at<float>(0,0) = camera[3];
    t.at<float>(1,0) = camera[4];
    t.at<float>(2,0) = camera[5];
    float f = camera[6];
    float k1 = camera[7];
    float k2 = camera[8];
    cv::Mat rotm = cv::Mat_<float>(3,3);
    cv::Rodrigues(rvec, rotm);//rodrigues transformation
    cv::Mat Xw = cv::Mat_<float>(3,1);
    Xw.at<float>(0,0) = point_3d[0];
    Xw.at<float>(1,0) = point_3d[1];
    Xw.at<float>(2,0) = point_3d[2];
    cv::Mat Xc = rotm * Xw + t;
    Xc = -Xc/Xc.at<float>(2,0);

    float l1 = cv::norm(Xc, cv::NORM_L2);
    float l2 = l1*l1;
    float rho = 1.0 + k1 * l1 + k2 * l2;
    cv::Mat Xp = f * rho * Xc;
    float d1, d2;
    d1 = Xp.at<float>(0,0) - observed_2d[0];
    d2 = Xp.at<float>(1,0) - observed_2d[1];
    return sqrt(d1*d1 + d2*d2);
}
/*******************************************************************************************************************************/
/* camera          : camera parameters to be optimized n*9 dimension matrix                                                     *
*  point_3d        : 3D points parameters to be optimized, n*3 dimension matrix                                                 *
*  observed_2d_data: Observed 2D data with dimension of n*2                                                                     *
*  cameraIndex     : Camera index, from 0 ... to number of cameras -1, which == camera.rows()-1                                 *
*  return value    : n * (9*nCam + 3 * nPoints) Jacobian matrix where nCam == number of camera, nPoints == number of 3D points  *
*********************************************************************************************************************************/
// cv::Mat Jacobian(Function, cv::Mat &cameras, cv::Mat & point_3ds,
//                             std::vector<unsigned> & cameraIndex, std::vector<unsigned> & pointIndex, 
//                             cv::Mat & observed_2d) {
//     //TODO ... The size of Jacobian Matrix should be nItem * (nCamera*9 + nPoints*3)
//     cv::Mat Jac = cv::Mat_<double>(pointIndex.size(), cameras.rows*9 + point_3ds.rows*3);
//     std::vector <double> cam;
//     std::vector<double> pt_3d;
//     std::vector<double> obs_2d;
//     unsigned nItem = cameraIndex.size();
//     for(unsigned i = 0; i< nItem; i++) {
//         cam.clear();
//         pt_3d.clear();
//         obs_2d.clear();
//         for(int j=0; j<N_CAMERA_PARAMETER; j++) {
//             cam.push_back(cameras.at<double>(cameraIndex[i], j));
//         }
//         //get point 3d parameters at index i
//         for(int j=0; j<N_POINT_3D; j++) {
//             pt_3d.push_back(point_3ds.at<double>(pointIndex[i], j));
//         }
//         for(int j=0; j<N_OBSERVED_2D; j++) {
//             obs_2d.push_back(observed_2d.at<double>(i,j));
//         }
//         //calculate jacobian matrix
         
//     }
//     //printf("Jacobian Matrix dimension = (%d, %d)\n",Jac.rows, Jac.cols);
//     return Jac;
// }

// calculate Jacobian vector for ith item
cv::Mat Jacobiani(Function func, cv::Mat & cameras, cv::Mat & point_3ds, 
                            std::vector<unsigned> & cameraIndex, std::vector<unsigned> &pointIndex,
                            cv::Mat & observed_2d , int item) {
    unsigned camera_id = cameraIndex[item];
    unsigned point_id = pointIndex[item];
    //get camera parameters and point parameters
    std::vector<float> camera_parameters;
    std::vector<float> point_parameters;
    std::vector<float> observation;
    camera_parameters.clear();
    point_parameters.clear();
    observation.clear();
    for(int i=0; i<9; i++) {
        camera_parameters.push_back(cameras.at<float>(camera_id, i));
    }
    for(int i=0; i<3; i++) {
        point_parameters.push_back(point_3ds.at<float>(point_id, i));
    }
    observation.push_back(observed_2d.at<float>(item,0));
    observation.push_back(observed_2d.at<float>(item,1));
    // item标明了Jacobian矩阵的行， camera_id*9+[0-8]标明了camera 参数的列， cameras.rows×9+point_id*3 + [0-2]表明了point的列， point在camera后面
    int column_dim = cameras.rows*9 + point_3ds.rows * 3;
    //printf("camera number = %i\n", column_dim);
    cv::Mat Jaci = cv::Mat_<float>(1, column_dim);
    for(int j=0; j < column_dim; j++)  {
        Jaci.at<float>(0,j) = 0.0;
    }
    //计算相机参数的第item行对应的 jacobian
    std::vector<float> cam1,cam2;
    cam1.clear();
    cam2.clear();
    for(int i=0; i<9; i++) {
        cam1.push_back(camera_parameters[i]);
        cam2.push_back(camera_parameters[i]);
    }
    for(int i=0; i<9; i++) {
        cam1[i] = cam1[i]+INFINITESIMAL;
        cam2[i] = cam1[i]-INFINITESIMAL;
        float jac = (func(cam1, point_parameters, observation) - func(cam2, point_parameters, observation))/(2.0 * INFINITESIMAL);
        Jaci.at<float>(0, camera_id*9 + i) = jac;
        cam1[i] = camera_parameters[i];
        cam2[i] = camera_parameters[i];
    }
    //计算3d点参数的第item行对应的3d point jacobian
    std::vector<float> pt1,pt2;
    pt1.clear();
    pt2.clear();
    for(int i=0; i < 3; i++) {
        pt1.push_back(point_parameters[i]);
        pt2.push_back(point_parameters[i]);
    }
    //printf("point parameters: \n");
    for(int i=0; i<3; i++) {
        //printf("debug info\n");
        //printf("pt1.size() = %lu\n", pt1.size());
        pt1[i] = pt1[i] + INFINITESIMAL;
        pt2[i] = pt2[i] - INFINITESIMAL;
        //printf("debug info\n");
        float jac = (func(camera_parameters, pt1, observation) - func(camera_parameters, pt2, observation))/(2.0 * INFINITESIMAL);
        Jaci.at<float>(0,cameras.rows*9 + point_id * 3 + i) = jac;
        //printf("%lf\n", )
        pt1[i] = point_parameters[i];
        pt2[i] = point_parameters[i];
    }
    return Jaci;
}
double TotalCost(Function, cv::Mat &cameras, cv::Mat & point_3ds,
                            std::vector<unsigned> & cameraIndex, std::vector<unsigned> & pointIndex, 
                            cv::Mat & observed_2d) {
    unsigned nItem = cameraIndex.size();
    
    double sum = 0.0;
    std::vector<float> cam;
    std::vector<float> pt_3d;
    std::vector<float> obs_2d;
    for(unsigned i = 0; i< nItem; i++ ) {
        cam.clear();
        pt_3d.clear();
        obs_2d.clear();
        //get camera parameters for camera[cam_idx]
        for(int j=0; j<N_CAMERA_PARAMETER; j++) {
            cam.push_back(cameras.at<float>(cameraIndex[i], j));
        }
        //get point 3d parameters at index i
        for(int j=0; j<N_POINT_3D; j++) {
            pt_3d.push_back(point_3ds.at<float>(pointIndex[i], j));
        }
        for(int j=0; j<N_OBSERVED_2D; j++) {
            obs_2d.push_back(observed_2d.at<float>(i,j));
        }
        float single_error = reprojection_err(cam, pt_3d, obs_2d);
        //printf("single error =%lf\n", single_error);
        sum += single_error;
        //TODO Jacobian i

    }
}
/*x[i+1] = x[i] - (H + s*I)^(-1)*f'(x[i])
 *so, J and f'(x) are needed for each step
 *question: how to select s. (turst region? line search?) 
 * refer to numerical optimization.pdf page 27, we can find the solution.
*/
cv::Mat Jacobian(Function func, cv::Mat & cameras, cv::Mat & point_3ds, 
                 std::vector<unsigned> & cameraIndex, std::vector<unsigned> &pointIndex,
                 cv::Mat & observed_2d) {
    unsigned nitem = cameraIndex.size();
    unsigned col = cameras.rows*9 + point_3ds.rows*3;
    cv::Mat Jac = cv::Mat::zeros(nitem, col, CV_32F);   
    for(int i=0; i<nitem; i++) {
        cv::Mat jaci = Jacobiani(func, cameras, point_3ds, cameraIndex, pointIndex, observed_2d, i);
        //printf("jaci.cols = %d\n", jaci.cols);
        for(int j = 0; j < jaci.cols; j++) {
            Jac.at<float>(i,j) = jaci.at<float>(0, j);
        }
    }
    return Jac;

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
    observe = cv::Mat_<float>(nObserve, 2);//n*2 matrix, n is number of observations
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
        float x,y;
        fscanf(fp, "%f",&x);
        fscanf(fp, "%f",&y);
        observe.at<float>(i, 0) = x;
        observe.at<float>(i, 1) = y;
    }
    //read the third part , camera parameters
    camera = cv::Mat_<float>(nCamera, 9);//n*9 dimension
    for(int i=0; i<nCamera; i++) {
        for(int j=0; j<9; j++) {
            float temp = 0.0;
            fscanf(fp, "%f", &temp);
            camera.at<float>(i,j) = temp;
        }
    }
    //read the forth part, point x,y,z
    points = cv::Mat_<float>(nPoint, 3);//n*3
    for(int i=0; i<nPoint; i++ ) {
        for(int j=0; j<3; j++) {
            float temp = 0.0;
            fscanf(fp, "%f", &temp);
            points.at<float>(i,j) = temp;
        }
    }
    fclose(fp);
}


#endif
