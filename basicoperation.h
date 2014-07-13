#ifndef __BASICOPERATION_H__
#define __BASICOPERATION_H__

#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <math.h>


using namespace cv;
using namespace std;

#define PI 3.141592653 
#define EPS 2.2204e-016 

Mat complexMatrixMul(Mat &m, Mat &n);
Mat complexMatrixDivide(Mat &m, Mat &n);

void imageInit(Mat &input, Mat &output);
void fft2(Mat &srcImage, Mat &dstImage);
void psf2otf(Mat &psf, Mat &otf, CvSize size);
void fftshift(Mat &input, Mat &output);
void fftShow(Mat &input);
void genaratePsf(Mat &psf, double len, double angle);
void conj(Mat &input, Mat &output);
void calConv(Mat &input, Mat &kernal, Mat &output);
void getpart(Mat &input, Mat &output, int i = 0);
void meshgrid(int x0, int xt, int y0, int yt, Mat &X, Mat &Y);
void calMag(Mat &input, Mat &output);

double estimatek(Mat &input);


#endif
