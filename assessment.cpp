#include "assessment.h"

/*****************************************************************
*函数名：calSpatialFreq
*功能：  计算图像的空间频率
*****************************************************************/
double  calSpatialFreq(Mat &input)
{
	double RF = 0, CF = 0, SF = 0;
	int i, j;
	int row = input.rows;
	int col = input.cols;
	int num = (row/-1)*(col-1);
	double s = 0;

	for (j = 0; j < row-1 ; j++)
	{
		uchar* current = input.ptr<uchar>(j);
		for (i = 0; i < col-1; i++)
		{
			s += (current[i + 1] - current[i])*(current[i + 1] - current[i]);
		}
	}
	RF = sqrt(s / num);

	for (j = 0; j < row-1 ; j++)
	{
		uchar* current = input.ptr<uchar>(j);
		uchar* next = input.ptr<uchar>(j + 1);
		for (i = 0; i < col-1; i++)
		{
			s += (next[i] - current[i])*(next[i] - current[i]);
		}
	}
	CF = sqrt(s / num);

	SF = sqrt(RF*RF + CF*CF);
	return SF;
}

/*****************************************************************
*函数名：calGMG
*功能：  计算图像灰度平均梯度值，其值越大表示图像越清晰
*****************************************************************/
double  calGMG(Mat &input)
{
	Mat I = Mat_<double>(input);
	int i, j;
	int row = input.rows;
	int col = input.cols;
	int num = (row -1)*(col - 1);
	double s = 0;
	double resolution = 0;

	for (j = 0; j < row-1; j++)
	{
		uchar* current = input.ptr<uchar>(j);
		uchar* next = input.ptr<uchar>(j + 1);
		for (i = 0; i < col-1; i++)
		{
			s += sqrt(((current[i + 1] - current[i])*(current[i + 1] - current[i]) + (next[i] - current[i])*(next[i] - current[i])) / 2);
		}
	}
	resolution = s / num;
	
	return resolution;
}


/********************************************************************************
*函数描述：	SpaceFreq 计算并返回一幅图像的空间频率
*函数参数：	IplImage *img 单通道8位图像
*函数返回值：double
*********************************************************************************/
double SpaceFreq(IplImage *img)
{
	double RF = 0;
	double CF = 0;
	double SF = 0;

	int i, j;//循环变量
	int height = img->height;
	int width = img->width;
	int step = img->widthStep / sizeof(uchar);
	uchar *data = (uchar*)img->imageData;
	double num = width*height;

	//行频计算
	for (i = 0; i<height; i++)
	{
		for (j = 0; j<width; j++)
		{
			RF += (data[i*step + j + 1] - data[i*step + j])*(data[i*step + j + 1] - data[i*step + j]);
		}
	}
	RF = sqrt(1.0*RF / num);

	//列频计算
	for (i = 0; i<height; i++)
	{
		for (j = 0; j<width; j++)
		{
			CF += (data[(i + 1)*step + j] - data[i*step + j])*(data[(i + 1)*step + j] - data[i*step + j]);
		}
	}
	CF = sqrt(1.0*CF / num);

	//空间频率
	SF = sqrt(RF*RF + CF*CF);
	return SF;
}
/********************************************************************************
*函数描述：	DefRto 计算并返回一幅图像的清晰度
*函数参数：	IplImage *img 单通道8位图像
*函数返回值：double
*********************************************************************************/
double DefRto(IplImage *img)
{
	double temp = 0;
	double DR = 0;
	int i, j;//循环变量
	int height = img->height;
	int width = img->width;
	int step = img->widthStep / sizeof(uchar);
	uchar *data = (uchar*)img->imageData;
	double num = width*height;

	for (i = 0; i<height; i++)
	{
		for (j = 0; j<width; j++)
		{
			temp += sqrt((pow((double)(data[(i + 1)*step + j] - data[i*step + j]), 2) + pow((double)(data[i*step + j + 1] - data[i*step + j]), 2)) / 2);
		}
	}
	DR = temp / num;
	return DR;
}


/*************************************************************************
*
* @函数名称：
*	calMSE()
*
* @输入参数:
*   Mat &image1           - 输入图像1
*	Mat &image2           - 输入图像2
*
* @返回值:
*   float                 - 返回计算出的MSE
*
* @输出：
*	无
*
* @说明:
*   该函数用来计算两幅图像的均方误差MSE
*
************************************************************************/
double calMSE(Mat &image1, Mat &image2)
{
	double MSE = 0, s = 0;
	int row1 = image1.rows;
	int col1 = image1.cols;
	int row2 = image2.rows;
	int col2 = image2.cols;
	int i = 0, j = 0;
	if (row1!=row2||col1!=col2)
		cerr<< "Image size not match!"<< endl;
	for (j = 0; j < row1; j++)
	{
		uchar* pvalue1 = image1.ptr<uchar>(j);
		uchar* pvalue2 = image2.ptr<uchar>(j);
		for (i = 0; i < col1; i++)
		{
			s+=(pvalue1[i] - pvalue2[i])*(pvalue1[i] - pvalue2[i]);
		}
	}
	MSE = s / (row1*col1);
	return MSE;
}

/*************************************************************************
*
* @函数名称：
*	calPSNR()
*
* @输入参数:
*   Mat &image1           - 输入图像1
*	Mat &image2           - 输入图像2
*
* @返回值:
*   float                 - 返回计算出的PSNR
*
* @输出：
*	无
*
* @说明:
*   该函数用来计算两幅图像的峰值信噪比PSNR，PSNR=10*log10((2^n-1)^2/MSE)
*
************************************************************************/
double calPSNR(Mat &image1, Mat &image2)
{
	double PSNR = 0, MSE = 0;

	MSE = calMSE(image1, image2);
	PSNR = 10 * log10(255 * 255 / MSE) ;

	return PSNR;
}

/*************************************************************************
*
* @函数名称：
*	calContrast()
*
* @输入参数:
*   Mat &image            - 输入图像
*
* @返回值:
*   float                 - 返回计算出的对比度
*
* @输出：
*	无
*
* @说明:
*   该函数用来计算单幅图像的对比度
*
************************************************************************/
double calContrast(Mat &image)
{
	double con = 0;
	double imin, imax;

	minMaxLoc(image, &imin, &imax);
	con = (imax - imin) / (imax + imin);


	return con;
}


/*************************************************************************
*
* @函数名称：
*	calWeight()
*
* @输入参数:
*   Mat &image            - 输入图像
*
* @返回值:
*   float                 - 返回计算出的对比度
*
* @输出：
*	无
*
* @说明:
*   该函数通过计算局部方差来获得各个像素点的权值，获得权值图像采用2X2块
*
************************************************************************/
void calWeight(Mat &input ,Mat &output)
{
	Mat temp = Mat_<double>(input);
	Mat mean ;

	int row = input.rows;
	int col = input.cols;

	for (int j = 0; j < row-2; j += 2)
	{
		for (int i = 0; i < col-2; i += 2)
		{
			Mat re(temp, Rect(i, j, 2, 2));
			meanStdDev(re,mean,re);
		}
	}
	temp = temp.mul(temp);
	Mat lc1(temp, Rect(col - 2, 0, 1, row));
	Mat lc2(temp, Rect(col - 1, 0, 1, row));

	Mat lr1(temp, Rect(0, row - 2, col, 1));
	Mat lr2(temp, Rect(0, row - 1, col, 1));

	lc2.copyTo(lc1);
	lr2.copyTo(lr1);

	output = temp.clone();

}
