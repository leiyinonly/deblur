#include "basicoperation.h"

/*************************************************************************
*
* @函数名称：
*	imageInit()
*
* @输入参数:
*   Mat &input             - 输入图像的引用
*
* @返回值:
*   无
*
* @输出：
*   Mat &input             - 输出图像的引用
*
* @说明:
*   该函数用来对图像尺寸进行预处理，对边界进行填充以便进行FFT
*
************************************************************************/
void imageInit(Mat &input,Mat &output)
{
	int m = getOptimalDFTSize(input.rows);
	int n = getOptimalDFTSize(input.cols);

	copyMakeBorder(input, output, 0, m - input.rows, 0, n - input.cols, BORDER_CONSTANT, Scalar::all(0));
}

/*************************************************************************
*
* @函数名称：
*	complexMatrixMul()
*
* @输入参数:
*   Mat &m             - 被乘数的引用，为一个两通道复数矩阵
*	Mat &n			   - 乘数的引用，为一个两通道复数矩阵
*
* @返回值:
*   Mat                - 返回复数乘法结果
*
* @输出：
*   无
*
* @说明:
*   该函数用来对复数进行点乘运算，输入和输出均为两通道矩阵，
*	两个通道分别存放实部和虚部
*
************************************************************************/
Mat complexMatrixMul(Mat &m, Mat &n)
{
	Mat a, b, c, d, result, resultRe, resultIm;
	Mat mp[] = { m, m };
	Mat np[] = { n, n };
	Mat resultp[] = { m, m };

	split(m, mp);
	split(n, np);
	a = mp[0];
	b = mp[1];
	c = np[0];
	d = np[1];

	resultRe = a.mul(c) - b.mul(d);
	resultIm = a.mul(d) + b.mul(c);

	resultp[0] = resultRe;
	resultp[1] = resultIm;

	merge(resultp, 2, result);

	return result;
}

/*************************************************************************
*
* @函数名称：
*	complexMatrixDivide()
*
* @输入参数:
*   Mat &m             - 被除数的引用，为一个两通道复数矩阵
*	Mat &n			   - 除数的引用，为一个两通道复数矩阵
*
* @返回值:
*   Mat                - 返回复数乘法结果
*
* @输出：
*   无
*
* @说明:
*   该函数用来对复数进行点除运算，输入和输出均为两通道矩阵，
*	两个通道分别存放实部和虚部
*
************************************************************************/
Mat complexMatrixDivide(Mat &m, Mat &n)
{
	Mat a, b, c, d, temp, result, resultRe, resultIm;
	Mat mp[] = { m, m };
	Mat np[] = { n, n };
	Mat resultp[] = { m, m };

	split(m, mp);
	split(n, np);
	a = mp[0];
	b = mp[1];
	c = np[0];
	d = np[1];

	temp = max(c.mul(c) + d.mul(d),EPS);
	divide((a.mul(c) + b.mul(d)), temp, resultRe);
	divide((b.mul(c) - a.mul(d)), temp, resultIm); 

	resultp[0] = resultRe;
	resultp[1] = resultIm;

	merge(resultp, 2, result);

	return result;
}

/*************************************************************************
*
* @函数名称：
*	psf2otf()
*
* @输入参数:
*   Mat &psf             - 待处理点扩散函数的引用
*   Mat &otf             - 用于存放输出的空间
*	CvSize size			 - 输出的尺寸
*
* @返回值:
*	无
*
* @输出：
*   Mat &otf             - 转换为指定大小的光学传递函数，为一个两通道矩阵
*
* @说明:
*   该函数用来把一个空间点扩散函数转换为指定大小的频域的光学传递函数
*
************************************************************************/
void psf2otf(Mat &psf, Mat &otf, CvSize size)
{
	Mat I = psf.clone();

	int m = getOptimalDFTSize(size.height);
	int n = getOptimalDFTSize(size.width);
	copyMakeBorder(psf, otf, 0, m - psf.rows, 0, n - psf.cols, BORDER_CONSTANT, Scalar::all(0));

	int cx = (I.cols-1) / 2;
	int cy = (I.rows-1) / 2;

	Mat q0(I, Rect(0, 0, cx, cy));   // Top-Left - Create a ROI per quadrant 
	Mat q1(I, Rect(cx, 0, cx, cy));  // Top-Right
	Mat q2(I, Rect(0, cy, cx, cy));  // Bottom-Left
	Mat q3(I, Rect(cx, cy, cx, cy)); // Bottom-Right

	Mat p0(otf, Rect(n-cx, m-cy, cx, cy));   // Top-Left - Create a ROI per quadrant 
	Mat p1(otf, Rect(0, m-cy, cx, cy));  // Top-Right
	Mat p2(otf, Rect(n-cx, 0, cx, cy));  // Bottom-Left
	Mat p3(otf, Rect(0, 0, cx, cy)); // Bottom-Right

	otf = Scalar::all(0);                           
	q0.copyTo(p0);
	q1.copyTo(p1);
	q2.copyTo(p2);
	q3.copyTo(p3);

	Mat temp = otf.clone();
	Mat planesOtf[] = { Mat_<float>(temp), Mat::zeros(temp.size(), CV_32F) };
	merge(planesOtf, 2, temp);
	dft(temp, otf);
	Mat t(otf, Rect(0, 0, size.width, size.height));
	t.copyTo(otf);
	//namedWindow("幅频谱图", 1);
	//imshow("幅频谱图", planesOtf[0]);
}

/******************************************************
*实现2维矩阵的离散傅里叶变换，返回一个两通道矩阵
******************************************************/
void fft2(Mat &srcImage, Mat &dstImage)
{
	Mat padded;
	//int m = getOptimalDFTSize(srcImage.rows);
	//int n = getOptimalDFTSize(srcImage.cols);

	//copyMakeBorder(srcImage, padded, 0, m - srcImage.rows, 0, n - srcImage.cols, BORDER_CONSTANT, Scalar::all(0));

	Mat planes[] = { Mat_<float>(srcImage), Mat::zeros(srcImage.size(), CV_32F) };

	merge(planes, 2, dstImage);
	dft(dstImage, dstImage);
}

/******************************************************
*将傅立叶变换后的数据交换象限便于显示
******************************************************/
void fftshift(Mat &input, Mat &output)
{
	Mat I = input(Rect(0, 0, input.cols & -2, input.rows & -2));

	int cx = I.cols / 2;
	int cy = I.rows / 2;

	Mat q0(I, Rect(0, 0, cx, cy));   // Top-Left - Create a ROI per quadrant 
	Mat q1(I, Rect(cx, 0, cx, cy));  // Top-Right
	Mat q2(I, Rect(0, cy, cx, cy));  // Bottom-Left
	Mat q3(I, Rect(cx, cy, cx, cy)); // Bottom-Right

	Mat tmp;                           // swap quadrants (Top-Left with Bottom-Right)
	q0.copyTo(tmp);
	q3.copyTo(q0);
	tmp.copyTo(q3);

	q1.copyTo(tmp);                    // swap quadrant (Top-Right with Bottom-Left)
	q2.copyTo(q1);
	tmp.copyTo(q2);

	output = I;
}

/******************************************************
*显示傅立叶变换后的对数幅度值
******************************************************/
void fftShow(Mat &input)
{
	Mat planes[] = { Mat::zeros(input.size(), CV_32F), Mat::zeros(input.size(), CV_32F) };
	split(input, planes);                   // planes[0] = Re(DFT(I), planes[1] = Im(DFT(I))
	magnitude(planes[0], planes[1], planes[0]);// planes[0] = magnitude  
	Mat magI = planes[0];

	magI += Scalar::all(1);                    // switch to logarithmic scale
	log(magI, magI);
	normalize(magI, magI, 1, 0, NORM_MINMAX);

	namedWindow("幅频谱图", 1);
	imshow("幅频谱图", magI);
}


/******************************************************
*计算复共轭
******************************************************/
void conj(Mat &input, Mat &output)
{
	Mat a, b;
	Mat mp[] = { input, input };
	Mat np[] = { input, input };

	split(input, mp);
	a = mp[0];
	b = mp[1];

	np[0] = a;
	np[1] = -b;

	merge(np, 2, output);
}

/******************************************************
*
*估计参数K，先通过计算局部方差估计信噪比
*局部方差的最大值为信号方差，最小值为噪声方差，再用经验公式修正
******************************************************/
double estimatek(Mat &input)
{
	double snr, k,max,min;
	CvSize size = input.size();
	Mat f = Mat_<float>(input),av,v;
	Mat g = f;// Mat(f, Rect(2, 2, size.width - 4, size.height - 4));
	Mat mask = Mat::ones(5, 5,CV_32F);
	mask = mask / 25;

	filter2D(g, av, -1, mask);
	pow((g - av), 2, v);
	filter2D(v, v, -1, mask);
	minMaxLoc(v,&min,&max);
	snr = 10 * log(max / min);
	//snr = 1.04*snr - 7;
	//k=pow(10, (-snr / 10));
	//k = 5 * k;

	return snr;
}

/******************************************************
*函数名：calConv
*功能：  用FFT计算卷积,时域卷积即频域相乘
******************************************************/
void calConv(Mat &input,Mat &kernal,Mat &output)
{
	CvSize size = input.size();
	Mat I, H,temp;
	fft2(input, I);
	psf2otf(kernal, H,size);
	temp = complexMatrixMul(I,H);
	dft(temp, output, DFT_INVERSE + DFT_SCALE);
}

/******************************************************
*函数名：getpart
*功能：  获取复数的实部或虚部，或者获取两通道矩阵的分量
******************************************************/
void getpart(Mat &input,  Mat &output,int i)
{
	Mat I = input.clone();
	Mat plane[] = { Mat_<float>(), Mat_<float>() };
	split(I, plane);
	if (i == 0)
		output = plane[0];
	else if (i==1)
		output = plane[1];
}

/*************************************************************************
*
* @函数名称：
*	genaratePsf()
*
* @输入参数:
*   Mat &psf               - 生成PSF的存储空间
*	double len             - 指定PSF的模糊长度
*	double angle           - 指定PSF的模糊角度
*
* @返回值:
*   无
*
* @输出：
*   Mat &psf               - 输出生成的模糊核
*
* @说明:
*   该函数用来生成指定参数的运动模糊核，模拟Matlab中的fspecial的‘motion’功能
*
************************************************************************/
void genaratePsf(Mat &psf, double len, double angle)
{
	len = max(1.0, len);
	double half = (len-1) / 2;
	double phi = angle/ 180.0* PI;
	double cosphi = cos(phi);
	double sinphi = sin(phi);
	double xsign;

	if (cosphi < 0)
	{
		xsign = -1;
	}
	else
	{
		if (cosphi == 90)
		{
			xsign = 0;
		}
		else
		{
			xsign = 1;
		}
	}

	double linewdt = 1.0;
	//计算0到90度PSF，实际角度-90到0度
	double sx = half*cosphi + linewdt*xsign - len*EPS;
	sx = cvFloor(sx);
	double sy = half*sinphi + linewdt - len*EPS;
	sy = cvFloor(sy);
	Mat X, Y;
	meshgrid(0, int(sx), 0, int(sy), X, Y);
	int row = X.rows;
	int col = X.cols;

	Mat_<float> dist2line=Mat::zeros(row, col, CV_32FC1);
	Mat_<float> rad(row, col, CV_32FC1);
	Mat_<float> h=Mat::zeros(row, col, CV_32FC1);

	addWeighted(Y, cosphi, X, -sinphi,0,dist2line);
	rad = X.mul(X) + Y.mul(Y);
	sqrt(rad,rad);

	/*为减小运算量，先计算一半大小的PSF*/
	for (int j = 0; j < row; j++)
	{
		float* prvalue = rad.ptr<float>(j);
		float* pdvalue = dist2line.ptr<float>(j);
		float* pxvalue = X.ptr<float>(j);
		float* pyvalue = Y.ptr<float>(j);

		for (int i = 0; i < col; i++)
		{
			if (prvalue[i] >= half && fabs(pdvalue[i]) <= linewdt)
			{
				uchar x = X.at<uchar>(j,i);
				float temp = float(half - fabs((x + pdvalue[i] * sinphi) / cosphi));
				pdvalue[i] = sqrt(pdvalue[i] * pdvalue[i] + temp*temp);
			}
			pdvalue[i] = float(linewdt + EPS - fabs(pdvalue[i]));
			if (pdvalue[i] < 0)
			{
				pdvalue[i] = 0;
			}
		}
	}

	/*将模糊核矩阵扩展至实际大小*/
	flip(dist2line, h, -1);
	int nrow = row + row - 1;
	int ncol = col + col - 1;
	Mat p = Mat::zeros(nrow, ncol, CV_32FC1);
	for (int j = 0; j < row; j++)
	{
		float* phvalue =h.ptr<float>(j);
		float* ppvalue = p.ptr<float>(j);
		for (int i = 0; i < col; i++)
		{
			ppvalue[i] = phvalue[i];
		}
	}

	for (int j = 0; j < row; j++)
	{
		float* ppvalue = p.ptr<float>(row-1+j);
		float* pdvalue = dist2line.ptr<float>(j);
		for (int i = 0; i < col; i++)
		{
			ppvalue[col-1+i] = pdvalue[i];
		}
	}

	/*保持图像总能量不变，归一化矩阵*/
	Scalar s1 = sum(p);
	p = p / s1[0];
	if (cosphi > 0)
	{
		flip(p, p, 0);
	}

	psf = p;
}

/*************************************************************************
*
* @函数名称：
*	meshgrid()
*
* @输入参数:
*   int x0           - 输入X轴的起始值
*	int xt           - 输入X轴的结束值
*	int y0           - 输入y轴的起始值
*	int yt           - 输入y轴的起始值
*	Mat &X           - 输出X的存储空间
*	Mat &Y           - 输出Y的存储空间
*
* @返回值:
*   无
*
* @输出：
*	Mat &X           - 生成的X矩阵
*	Mat &Y           - 生成的Y矩阵
*
* @说明:
*   该函数用来模拟Matlab中的meshgrid函数
*
************************************************************************/
void meshgrid(int x0, int xt, int y0, int yt, Mat &X, Mat &Y)
{
	int i, j;
	int jd = yt + 1 - y0;
	int id = xt + 1 - x0;
	X = Mat(jd, id, CV_8UC1);
	Y = Mat(jd, id, CV_8UC1);
	
	for (j = 0; j < jd; j++)
	{
		uchar* pxvalue = X.ptr<uchar>(j);
		uchar* pyvalue = Y.ptr<uchar>(j);
		for (i = 0; i < id; i++)
		{
			pxvalue[i] = x0 + i;
			pyvalue[i] = y0 + j;
		}
	}

}

//计算两通道复数的幅值
void calMag(Mat &input, Mat &output)
{
	Mat planes[] = { Mat::zeros(input.size(), CV_32F), Mat::zeros(input.size(), CV_32F) };

	split(input, planes);
	magnitude(planes[0], planes[1], planes[0]);
	output = planes[0];
}
