#ifndef __RESTORE_H__
#define __RESTORE_H__

#include"basicoperation.h"

void inverseFilter(Mat &input, Mat &kernal, Mat &output);
void wnrFilter1(Mat &input, Mat &kernal, Mat &output);
void wnrFilter2(Mat &input, Mat &kernal, Mat &output);
void deconvRL1(Mat &input, Mat &kernal, Mat &output, int num=10);
void deconvRL2(Mat &input, Mat &kernal, Mat &output, int num = 10);
void deconvRL3(Mat &input, Mat &kernal, Mat &output, int num = 10);
void deconvRL4(Mat &input, Mat &kernal, Mat &output, int num = 10);
void deconvRL5(Mat &input, Mat &kernal, Mat &output, int num);

Mat corelucy(Mat &Y, Mat &H, Mat &g);
Mat edgetaper(Mat &image, Mat &psf);

#endif

