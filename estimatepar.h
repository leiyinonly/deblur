#ifndef __ESTIMATEPAR_H__
#define __ESTIMATEPAR_H__

#include"basicoperation.h"

double  directionDiffAngle(Mat &input);
double  directionDiffAngle1(Mat &input);

void sectionMask(double angle, Mat &mask);
void sectionMask1(double angle, Mat &mask);
void cepstral(Mat &input, Mat &output);
//void matrix_orient(float angle, double *motion_matrix);
//char get_angle(double *motion_matrix, unsigned char*data,
//	const unsigned short width, const unsigned short height);
//unsigned char motion_orient(unsigned char *Y, unsigned int i, unsigned int j,
//double *motion_matrix, float angle, unsigned int width);



#endif

