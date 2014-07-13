#include "estimatepar.h"

/*****************************************************************
*函数名：directionDiffAngle
*功能：  方向微分法估计运动模糊角度,微元半径为2，3X3微分算子
*****************************************************************/
double  directionDiffAngle(Mat &input)
{
	//3X3算子
	Mat mask = Mat_<double>(3, 3);
	double angle=0,step=1;
	float n=0;

	double s0 = 255;
	double sx = 0,sy=0,temp;
	int i, j;
	int row = input.rows;
	int col = input.cols;

	for (angle = -90; angle <= 90; angle = angle + step)
	{
		sx = 0;
		sectionMask(angle, mask);

		Scalar ss = sum(mask);

		for (j = 2; j < row - 2; j=j+4)
		{
			sy = 0;
			uchar* previous1 = input.ptr<uchar>(j - 1);
			uchar* previous2 = input.ptr<uchar>(j - 2);
			uchar* current = input.ptr<uchar>(j);
			uchar* next1 = input.ptr<uchar>(j + 1);
			uchar* next2 = input.ptr<uchar>(j + 2);
			for (i = 0; i < col - 2; i = i + 2)
			{
				if (angle >= -90.0&&angle < 0.0)
				{
					temp =(mask.at<double>(0)*previous2[i] + mask.at<double>(1)*previous2[i + 1] + mask.at<double>(2)*previous2[i + 2] +
						mask.at<double>(3)*previous1[i] + mask.at<double>(4)*previous1[i + 1] + mask.at<double>(5)*previous1[i + 2] +
						mask.at<double>(6)*current[i] + mask.at<double>(7)*current[i + 1] + mask.at<double>(8)*current[i + 2]);
					if (temp>255)
						temp = 255;
					temp = temp / (col - 2);
				}
				else
				{
					temp = (mask.at<double>(0)*current[i] + mask.at<double>(1)*current[i + 1] + mask.at<double>(2)*current[i + 2] +
						mask.at<double>(3)*next1[i] + mask.at<double>(4)*next1[i + 1] + mask.at<double>(5)*next1[i + 2] +
						mask.at<double>(6)*next2[i] + mask.at<double>(7)*next2[i + 1] + mask.at<double>(8)*next2[i + 2]);
					if (temp>255)
						temp = 255;
					temp = temp / (col - 2);
				}
				sy = sy + temp;
			}
			sx += sy / (row - 4);
		}

		if (sx < s0)
		{
			s0 = sx;
			n = angle;
		}
	}
	angle = -angle;
	if (angle < 0)
		angle = 180+angle;

	return angle;
}


/*****************************************************************
*函数名：directionDiffAngle1
*功能：  方向微分法估计运动模糊角度,微元半径为1，2X2微分算子
*****************************************************************/
double  directionDiffAngle1(Mat &input)
{
	Mat mask = Mat_<double>(2, 2);
	double angle = 0, step = 1;
	int m = 0, n = 0;

	double s0 = 255;
	double sx = 0, sy = 0, temp;
	int i, j;
	int row = input.rows;
	int col = input.cols;


	for (angle = -90; angle <= 90; angle = angle + step)
	{
		sx = 0;
		sectionMask1(angle, mask);
		for (j = 1; j < row - 1; j = j + 1)
		{
			sy = 0;
			uchar* previous = input.ptr<uchar>(j - 1);
			uchar* current = input.ptr<uchar>(j);
			uchar* next = input.ptr<uchar>(j + 1);
			for (i = 0; i < col - 1; i = i + 1)
			{
				if (angle >= -90.0&&angle < 0.0)
				{
					temp = (mask.at<double>(0)*previous[i] + mask.at<double>(1)*previous[i + 1] +
							mask.at<double>(2)*current[i] + mask.at<double>(3)*current[i + 1]) ;
					if (temp>255)
						temp = 255;
					temp = temp / (col - 1);
				}
				else
				{
					temp = (mask.at<double>(0)*current[i] + mask.at<double>(1)*current[i + 1] +
						mask.at<double>(2)*next[i] + mask.at<double>(3)*next[i + 1]);
					if (temp>255)
						temp = 255;
					temp = temp / (col - 1);
				}
				sy = sy + temp;
			}
			sx += sy / (row - 2);
		}

		if (sx < s0)
		{
			s0 = sx;
			n = m;
		}
		m++;
	}
	angle = -90 + n*step;
	if (angle > 0)
		angle = 180 - angle;
	angle = -angle;

	return angle;
}
/*****************************************************************
*函数名：sectionMask1
*功能：  判定待插值像素区域，计算双线性插值后方向微分掩膜
*****************************************************************/
void sectionMask1(double angle, Mat &mask)
{
	if (angle >= -90.0&&angle < 0)
	{
		angle = angle*PI / 180.0;
		mask.at<double>(0) = 1-cos(angle) + sin(angle) - sin(angle)*cos(angle);
		mask.at<double>(1) = cos(angle) + sin(angle)*cos(angle);
		mask.at<double>(2) = -sin(angle) - sin(angle)*cos(angle)-1;
		mask.at<double>(3) = -sin(angle)*cos(angle);
	}
	else if (angle >= 0 && angle < 90.0)
	{
		angle = angle*PI / 180.0;
		mask.at<double>(0) = -cos(angle) + sin(angle) - sin(angle)*cos(angle);
		mask.at<double>(1) = cos(angle) + sin(angle)*cos(angle);
		mask.at<double>(2) = -sin(angle)  - sin(angle)*cos(angle);
		mask.at<double>(3) = - sin(angle)*cos(angle);
	}
}

/*****************************************************************
*函数名：sectionMask
*功能：  判定待插值像素区域，计算双线性插值后方向微分掩膜
*****************************************************************/
void sectionMask(double angle, Mat &mask)
{
	if (angle >= -90.0&&angle<-60.0)
	{
		angle = angle*PI / 180.0;
		mask.at<double>(0) = -1 - 2 * sin(angle) + 2 * cos(angle) + 4 * sin(angle)*cos(angle);
		mask.at<double>(1) = -2 * cos(angle) - 4 * sin(angle)*cos(angle);
		mask.at<double>(2) = 0.0;
		mask.at<double>(3) = 2 + 2 * sin(angle) - 4 * cos(angle) - 4 * sin(angle)*cos(angle);
		mask.at<double>(4) = 4 * cos(angle) + 4 * sin(angle)*cos(angle);
		mask.at<double>(5) = 0.0;
		mask.at<double>(6) = -1.0;
		mask.at<double>(7) = 0.0;
		mask.at<double>(8) = 0.0;
	}
	else if (angle >= -60.0&&angle<-30.0)
	{
		angle = angle*PI / 180.0;
		mask.at<double>(0) = 0.0;
		mask.at<double>(1) = -2 - 4 * sin(angle) + 2 * cos(angle) + 4 * sin(angle)*cos(angle);
		mask.at<double>(2) = 1 + 2 * sin(angle) - 2 * cos(angle) - 4 * sin(angle)*cos(angle);
		mask.at<double>(3) = 0.0;
		mask.at<double>(4) = 4 + 4 * sin(angle) - 4 * cos(angle) - 4 * sin(angle)*cos(angle);
		mask.at<double>(5) = -2 - 2 * sin(angle) + 4 * cos(angle) + 4 * sin(angle)*cos(angle);
		mask.at<double>(6) = -1.0;
		mask.at<double>(7) = 0.0;
		mask.at<double>(8) = 0.0;
	} 
	else if (angle >= -30.0&&angle<0.0)
	{
		angle = angle*PI / 180.0;
		mask.at<double>(0) = 0.0;
		mask.at<double>(1) = 0.0;
		mask.at<double>(2) = 0.0;
		mask.at<double>(3) = 0.0;
		mask.at<double>(4) = -4 * sin(angle) + 4 * sin(angle)*cos(angle);
		mask.at<double>(5) = 2 * sin(angle) - 4 * sin(angle)*cos(angle);
		mask.at<double>(6) = -1.0;
		mask.at<double>(7) = 2 + 4 * sin(angle) - 2 * cos(angle) - 4 * sin(angle)*cos(angle);
		mask.at<double>(8) = -1 - 2 * sin(angle) + 2 * cos(angle) + 4 * sin(angle)*cos(angle);
	}
	else if (angle >= 0.0&&angle<30.0)
	{
		angle = angle*PI / 180.0;
		mask.at<double>(0) = -1.0;
		mask.at<double>(1) = 2 - 4 * sin(angle) - 2 * cos(angle) + 4 * sin(angle)*cos(angle);
		mask.at<double>(2) = -1 + 2 * sin(angle) + 2 * cos(angle) - 4 * sin(angle)*cos(angle);
		mask.at<double>(3) = 0.0;
		mask.at<double>(4) = 4 * sin(angle) - 4 * sin(angle)*cos(angle);
		mask.at<double>(5) = -2 * sin(angle) + 4 * sin(angle)*cos(angle);
		mask.at<double>(6) = 0.0;
		mask.at<double>(7) = 0.0;
		mask.at<double>(8) = 0.0;
	}
	else if (angle >= 30.0&&angle<60.0)
	{
		angle = angle*PI / 180.0;
		mask.at<double>(0) = -1.0,
		mask.at<double>(1) = 0.0;
		mask.at<double>(2) = 0.0;
		mask.at<double>(3) = 0.0;
		mask.at<double>(4) = 4 - 4 * sin(angle) - 4 * cos(angle) + 4 * sin(angle)*cos(angle);
		mask.at<double>(5) = -2 + 2 * sin(angle) + 4 * cos(angle) - 4 * sin(angle)*cos(angle);
		mask.at<double>(6) = 0.0;
		mask.at<double>(7) = -2 + 4 * sin(angle) + 2 * cos(angle) - 4 * sin(angle)*cos(angle);
		mask.at<double>(8) = 1 - 2 * sin(angle) - 2 * cos(angle) + 4 * sin(angle)*cos(angle);
	}
	else
	{
		angle = angle*PI / 180.0;
		mask.at<double>(0) = -1.0;
		mask.at<double>(1) = 0.0;
		mask.at<double>(2) = 0.0;
		mask.at<double>(3) = 2 - 2 * sin(angle) - 4 * cos(angle) + 4 * sin(angle)*cos(angle);
		mask.at<double>(4) = 4 * cos(angle) - 4 * sin(angle)*cos(angle);
		mask.at<double>(5) = 0.0;
		mask.at<double>(6) = -1 + 2 * sin(angle) + 2 * cos(angle) - 4 * sin(angle)*cos(angle);
		mask.at<double>(7) = -2 * cos(angle) + 4 * sin(angle)*cos(angle);
		mask.at<double>(8) = 0.0;
	}
}


void cepstral(Mat &input, Mat &output)
{
	Mat lg, lgcep,temp;

	fft2(input, temp);
	calMag(temp, temp);
	log((1 + temp),lg );
	dft(lg, lgcep, DFT_INVERSE + DFT_SCALE);
	fftshift(lgcep, lgcep);
	calMag(lgcep, temp);
	output = temp;
}