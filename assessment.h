#ifndef __ASSESSMENT_H__
#define __ASSESSMENT_H__

#include"basicoperation.h"

double  calSpatialFreq(Mat &input);
double  calGMG(Mat &input);
double SpaceFreq(IplImage *img);
double DefRto(IplImage *img);
double calMSE(Mat &image1, Mat &image2);
double calPSNR(Mat &image1, Mat &image2);
double calContrast(Mat &image);
void calWeight(Mat &input, Mat &output);

#endif
