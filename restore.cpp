#include "restore.h"

/*****************************************************************
*函数名：inverseFilter
*功能：  实现逆滤波F(u,v)=G(u,v)/H(u,v)
*****************************************************************/
void inverseFilter(Mat &input,Mat &kernal, Mat &output)
{
	Mat F, H, G,temp;
	Mat plane[] = { Mat_<float>(input), Mat_<float>(input) };

	fft2(input, G);
	psf2otf(kernal, H,input.size());
	H = H + 0.05;
	F = complexMatrixDivide(G,H);
	dft(F, temp, DFT_INVERSE + DFT_SCALE);
	split(temp, plane);
	magnitude(plane[0], plane[1], plane[0]);// planes[0] = magnitude  
	plane[0] = Mat_<uchar>(plane[0]);
	output = plane[0].clone();
}

/*****************************************************************
*函数名：wnrFilter1
*功能：  实现维纳滤波
*F(u,v)=G(u,v)*(conj(H(u,v)/(|H(u,v)|*|H(u,v)|+a))
*F(u,v)=(G(u,v)/H(u,v))*|H(u,v)|*|H(u,v)|/(|H(u,v)|*|H(u,v)|+a)
*****************************************************************/
void wnrFilter1(Mat &input, Mat &kernal, Mat &output)
{
	Mat input1 = edgetaper(input, kernal);
	Mat F, H, HJ, HH, G;
	fft2(input1, G);
	CvSize size = input1.size();
	Mat plane[] = { Mat_<float>(input1), Mat::zeros(input1.size(), CV_32F) };
	psf2otf(kernal, H, size);
	conj(H, HJ);
	HH = complexMatrixMul(H, HJ);
//	double k = estimatek(input);
	HH = HH/(HH + 0.05);
	F = complexMatrixDivide(G, H);
	F = complexMatrixMul(HH, F);
	dft(F, F, DFT_INVERSE + DFT_SCALE);
	
	Mat FF, FJ;
	conj(F, FJ);
	FF = complexMatrixMul(F, FJ);
	split(FF, plane);
	FF = plane[0];
	split(F, plane);
	plane[0] = ((plane[0] + FF)/2);
	plane[0] = Mat_<uchar>(plane[0]);
	output = plane[0].clone();
}

/*****************************************************************
*函数名：wnrFilter2
*功能：  实现改进的维纳滤波,采用边界延拓抑制振铃
*F(u,v)=G(u,v)*(conj(H(u,v)/(|H(u,v)|*|H(u,v)|+a))
*F(u,v)=(G(u,v)/H(u,v))*|H(u,v)|*|H(u,v)|/(|H(u,v)|*|H(u,v)|+a)
*****************************************************************/
void wnrFilter2(Mat &input, Mat &kernal, Mat &output)
{
	CvSize size = input.size();
	CvSize kernalsize = kernal.size();
	CvSize nsize;
	
	nsize.width=2*size.width;
	nsize.height = 2 * size.height;
	
	Mat y = input.clone(),y1,y2;
	Mat yy = Mat::zeros(nsize, CV_8U);
	Mat yy0 = Mat(yy, Rect(0, 0, size.width, size.height));
	Mat yy1 = Mat(yy, Rect(size.width, 0, size.width, size.height));
	Mat yy2 = Mat(yy, Rect(0, 0, nsize.width, size.height));
	Mat yy3 = Mat(yy, Rect(0, size.height, nsize.width, size.height));
	
	y.copyTo(yy0);
	flip(y, y1, 1);
	y1.copyTo(yy1);
	flip(yy2, y2, 0);
	y2.copyTo(yy3);
	
	Mat F, H, HJ, HH, G;
	Mat plane[] = { Mat_<float>(yy), Mat::zeros(yy.size(), CV_32F) };
	
	fft2(yy, G);
	psf2otf(kernal, H, nsize);
	conj(H, HJ);
	HH = complexMatrixMul(H, HJ);
//	double k = estimatek(input);
	HH = HH/(HH + 0.05);
	F = complexMatrixDivide(G, H);
	F = complexMatrixMul(HH, F);
	dft(F, F, DFT_INVERSE + DFT_SCALE);
			
	Mat FF, FJ;
	conj(F, FJ);
	FF = complexMatrixMul(F, FJ);
	split(FF, plane);
	FF = plane[0];
	split(F, plane);
	plane[0] = ((plane[0] + FF)/2);
	plane[0] = Mat(plane[0], Rect(0, 0, size.width, size.height));
	plane[0] = Mat_<uchar>(plane[0]);
	output = plane[0].clone();
}


/*****************************************************************
*函数名：deconvRL1
*功能：  实现RL迭代复原,使用傅里叶变换
*f(k+1)=f(k)(corr((g(k)/(conv(h,f(k)))),hInv)
*****************************************************************/
void deconvRL1(Mat &input, Mat &kernal, Mat &output, int num)
{
	Mat J = Mat_<float>(input);
	Mat psfInv, g;
	CvSize size = input.size();

	g = Mat::ones(size, CV_32F);
	g = 0.5*g;
	flip(kernal, psfInv, -1);

	Mat H,H1,G;
	psf2otf(kernal, H, size);
	psf2otf(psfInv, H1, size);

	//Scalar a = sum(psfInv);

	for (int i = 0; i < num; i++)
	{
		Mat RE,temp;
		calConv(g, kernal, temp);
		//fft2(g, G);
		//temp = complexMatrixMul(H, G);
		//dft(temp, temp, DFT_INVERSE + DFT_SCALE);
		//	cout << "min=" << temp << endl;
		getpart(temp, RE);
		max(RE, EPS, RE);
		divide(J, RE, temp);
		calConv(temp, psfInv, temp);
		//fft2(temp, temp);
		//temp = complexMatrixMul(H1, temp);
		//dft(temp, temp, DFT_INVERSE + DFT_SCALE);

		getpart(temp, RE);
		g = g.mul(RE);
	}
	//double min, max;
	//minMaxLoc(g, &min, &max);
	//cout <<"min="<< min<<"max="<< max<< endl;
	g = Mat_<uchar>(g);
	//normalize(g, g, 1, 0, NORM_MINMAX);

	output = g.clone();
}

/*****************************************************************
*函数名：deconvRL2
*功能：  实现RL迭代复原,使用傅里叶变换,进行边界延拓
*f(k+1)=f(k)(corr((g(k)/(conv(h,f(k)))),hInv)
*****************************************************************/
void deconvRL2(Mat &input, Mat &kernal, Mat &output, int num)
{
	CvSize size = input.size();
	CvSize kernalsize = kernal.size();
	CvSize nsize;

	nsize.width = 2 * size.width;
	nsize.height = 2 * size.height;

	Mat y = input.clone(), y1, y2;
	Mat yy = Mat::zeros(nsize, CV_8U);
	Mat yy0 = Mat(yy, Rect(0, 0, size.width, size.height));
	Mat yy1 = Mat(yy, Rect(size.width, 0, size.width, size.height));
	Mat yy2 = Mat(yy, Rect(0, 0, nsize.width, size.height));
	Mat yy3 = Mat(yy, Rect(0, size.height, nsize.width, size.height));

	y.copyTo(yy0);
	flip(y, y1, 1);
	y1.copyTo(yy1);
	flip(yy2, y2, 0);
	y2.copyTo(yy3);

	Mat J = Mat_<float>(yy);
	Mat psfInv, g;

	g = Mat::ones(nsize, CV_32F);
	g = 0.5*g;
	flip(kernal, psfInv, -1);

	Mat H, H1, G;
	psf2otf(kernal, H, nsize);
	psf2otf(psfInv, H1, nsize);

	//Scalar a = sum(psfInv);

	for (int i = 0; i < num; i++)
	{
		Mat RE, temp;
		calConv(g, kernal, temp);
		//fft2(g, G);
		//temp = complexMatrixMul(H, G);
		//dft(temp, temp, DFT_INVERSE + DFT_SCALE);
		//	cout << "min=" << temp << endl;
		getpart(temp, RE);

		max(RE, EPS, RE);
		divide(J, RE, temp);
		calConv(temp, psfInv, temp);
		//fft2(temp, temp);
		//temp = complexMatrixMul(H1, temp);
		//dft(temp, temp, DFT_INVERSE + DFT_SCALE);

		getpart(temp, RE);
		g = g.mul(RE);
		//	if (i == 3) 		cout << "min=" << RE << endl;
	}
	g = Mat(g, Rect(0, 0, size.width, size.height));

	//double min, max;
	//minMaxLoc(g, &min, &max);
	//cout << "min=" << min << "max=" << max << endl;
	g = Mat_<uchar>(g);

	output = g.clone();
}


/*****************************************************************
*函数名：deconvRL3
*功能：  实现RL迭代复原，使用FFT计算卷积，使用矢量外推法进行优化加速
*f(k+1)=f(k)(corr((g(k)/(conv(h,f(k)))),hInv)
*      =beta(f(k))
*y(k)=x(k)+lambda(k)*h(k)
*h(k)=x(k)-x(k-1)
*x(k+1)=y(k)+g(k)
*g(k)=beta(y(k))-y(k)
*lambda(k)=sum(g(k-1).*g(k-2))/sum(g(k-2).*g(k-2)),0<lambda(k)<1
*****************************************************************/
void deconvRL3(Mat &input, Mat &kernal, Mat &output,int num)
{
	Mat input1 = edgetaper(input, kernal);
	Mat psfInv, g;
	CvSize size = input.size();

	g = Mat::ones(size, CV_32F);
	g = 0.5*g;
	flip(kernal, psfInv, -1);

	Mat J[] = { Mat_<float>(input1), Mat_<float>(g), Mat::zeros(input1.size(), CV_32F) };
	Mat L[] = { Mat::zeros(input1.size(), CV_32F), Mat::zeros(input1.size(), CV_32F) };
	Mat H, H1, G, Y;
	float lambda = 0;
	psf2otf(kernal, H, size);
	psf2otf(psfInv, H1, size);
	Scalar s1, s2;

	for (int i = 0; i < num; i++)
	{
		Mat temp, g1, g2;
		Mat TEMP[] = { Mat::zeros(size, CV_32FC2), Mat::zeros(size, CV_32FC2) };

		if (i>1)
		{
			g1 = L[0].mul(L[1]);
			g2 = L[1].mul(L[1]);
			s1 = sum(g1);
			s2 = sum(g2);
			lambda = float(s1[0] / max(s2[0] ,EPS));
			lambda = max(min(lambda, float(1.0)), float(0));
		}

		Y = J[1] + lambda*(J[1] - J[2]);
		Y = max(Y,0);

		Mat CC = corelucy(Y, H, J[0]);

		temp = complexMatrixMul(CC, H1);
		dft(temp, temp, DFT_INVERSE + DFT_SCALE);
		split(temp, TEMP);

		J[2] = J[1].clone();
		J[1] = TEMP[0].mul(Y);

		L[1] = L[0].clone();
		L[0] = J[1] - Y;
	}
	J[1] = Mat_<uchar>(J[1]);                                       
	output = J[1].clone();
}

/*****************************************************************
*函数名：deconvRL4
*功能：  实现RL迭代复原，使用FFT计算卷积，使用矢量外推法进行优化加速
*f(k+1)=f(k)(corr((g(k)/(conv(h,f(k)))),hInv)
*      =beta(f(k))
*y(k)=x(k)+lambda(k)*h(k)
*h(k)=x(k)-x(k-1)
*x(k+1)=y(k)+g(k)
*g(k)=beta(y(k))-y(k)
*lambda(k)=sum(g(k-1).*g(k-2))/sum(g(k-2).*g(k-2)),0<lambda(k)<1
*****************************************************************/
void deconvRL4(Mat &input, Mat &kernal, Mat &output, int num)
{
	CvSize size = input.size();
	CvSize kernalsize = kernal.size();
	CvSize nsize;

	nsize.width = 2 * size.width;
	nsize.height = 2 * size.height;

	Mat y = input.clone(), y1, y2;
	Mat yy = Mat::zeros(nsize, CV_8U);
	Mat yy0 = Mat(yy, Rect(0, 0, size.width, size.height));
	Mat yy1 = Mat(yy, Rect(size.width, 0, size.width, size.height));
	Mat yy2 = Mat(yy, Rect(0, 0, nsize.width, size.height));
	Mat yy3 = Mat(yy, Rect(0, size.height, nsize.width, size.height));

	y.copyTo(yy0);
	flip(y, y1, 1);
	y1.copyTo(yy1);
	flip(yy2, y2, 0);
	y2.copyTo(yy3);

	Mat psfInv, g;

	g = Mat::ones(nsize, CV_32F);
	g = 0.5*g;
	flip(kernal, psfInv, -1);

	Mat J[] = { Mat_<float>(yy), Mat_<float>(g), Mat::zeros(nsize, CV_32F) };
	Mat L[] = { Mat::zeros(nsize, CV_32F), Mat::zeros(nsize, CV_32F) };
	Mat H, H1, G, Y;
	float lambda = 0;
	psf2otf(kernal, H, nsize);
	psf2otf(psfInv, H1, nsize);
	Scalar s1, s2;

	for (int i = 0; i < num; i++)
	{
		Mat temp, g1, g2;
		Mat TEMP[] = { Mat::zeros(size, CV_32FC2), Mat::zeros(size, CV_32FC2) };

		if (i>1)
		{
			g1 = L[0].mul(L[1]);
			g2 = L[1].mul(L[1]);
			s1 = sum(g1);
			s2 = sum(g2);
			lambda = float(s1[0] / max(s2[0], EPS));
			lambda = max(min(lambda, float(1.0)), float(0));
		}

		Y = J[1] + lambda*(J[1] - J[2]);
		Y = max(Y, 0);

		Mat CC = corelucy(Y, H, J[0]);

		temp = complexMatrixMul(CC, H1);
		dft(temp, temp, DFT_INVERSE + DFT_SCALE);
		split(temp, TEMP);

		J[2] = J[1].clone();
		J[1] = TEMP[0].mul(Y);

		L[1] = L[0].clone();
		L[0] = J[1] - Y;
	}
	J[1] = Mat(J[1], Rect(0, 0, size.width, size.height));
	J[1] = Mat_<uchar>(J[1]);
	output = J[1].clone();
}

/*****************************************************************
*函数名：deconvRL5
*功能：  实现RL迭代复原，使用FFT计算卷积，使用矢量外推法进行优化加速
*f(k+1)=f(k)(corr((g(k)/(conv(h,f(k)))),hInv)
*      =beta(f(k))
*y(k)=x(k)+lambda(k)*h(k)
*h(k)=x(k)-x(k-1)
*x(k+1)=y(k)+g(k)
*g(k)=beta(y(k))-y(k)
*lambda(k)=sum(g(k-1).*g(k-2))/sum(g(k-2).*g(k-2)),0<lambda(k)<1
*****************************************************************/
void deconvRL5(Mat &input, Mat &kernal, Mat &output, int num)
{
	Mat input1 = edgetaper(input, kernal);
	Mat psfInv, g;
	CvSize size = input.size();

	g = Mat::ones(size, CV_32F);
	g = 0.5*g;
	flip(kernal, psfInv, -1);

	Mat J[] = { Mat_<float>(input1), Mat_<float>(g), Mat::zeros(input1.size(), CV_32F) };
	Mat L[] = { Mat::zeros(input1.size(), CV_32F), Mat::zeros(input1.size(), CV_32F) };
	Mat H, H1, G, Y;
	float lambda = 0;
	psf2otf(kernal, H, size);
	psf2otf(psfInv, H1, size);
	Scalar s1, s2;

	for (int i = 0; i < num; i++)
	{
		Mat temp, g1, g2;
		Mat TEMP[] = { Mat::zeros(size, CV_32FC2), Mat::zeros(size, CV_32FC2) };

		if (i>1)
		{
			g1 = L[0].mul(L[1]);
			g2 = L[1].mul(L[1]);
			s1 = sum(g1);
			s2 = sum(g2);
			lambda = float(s1[0] / max(s2[0], EPS));
			lambda = max(min(lambda, float(1.0)), float(0));
		}

		Y = J[1] + lambda*(J[1] - J[2]);
		Y = max(Y, 0);

		Mat CC = corelucy(Y, H, J[0]);

		temp = complexMatrixMul(CC, H1);
		dft(temp, temp, DFT_INVERSE + DFT_SCALE);
		split(temp, TEMP);

		J[2] = J[1].clone();
		J[1] = TEMP[0].mul(Y);

		Mat gx =(Mat_<int>(1,2)<<-1, 1);
		filter2D(Y, Y, -1, gx);
		double l1=norm(Y, NORM_L1);
		double I = 1/(1 +0.002*0.002*0.002*l1/0.005);
		cout << "I=" << I<< endl;
		J[1] = J[1] *0.98;

		L[1] = L[0].clone();
		L[0] = J[1] - Y;
	}
	J[1] = Mat_<uchar>(J[1]);
	output = J[1].clone();
}

Mat corelucy(Mat &Y,Mat &H ,Mat &g)
{
	Mat YF, temp;
	Mat TEMP[] = { Mat::zeros(Y.size(), CV_32FC2), Mat::zeros(Y.size(), CV_32FC2) }; 

	fft2(Y, YF);
	temp = complexMatrixMul(H, YF);
	dft(temp, temp, DFT_INVERSE + DFT_SCALE);
	split(temp, TEMP);
	temp = g / TEMP[0];
	fft2(temp, temp);

	return temp;
}


Mat edgetaper(Mat &image, Mat &psf)
{
	int irow = image.rows;
	int icol = image.cols;
	int prow = psf.rows;
	int pcol = psf.cols;
	int temp=0;

	Mat blur;
	Mat psfRowProj(1, prow, CV_32FC1);
	Mat psfColProj(1, pcol, CV_32FC1);
	Mat planes1[] = { Mat::zeros(image.size(), CV_32FC1), Mat::zeros(image.size(), CV_32FC1) };
	int i, j;

	calConv(image, psf, blur);
	split(blur, planes1);
	blur = planes1[0];

	for (j = 0; j < prow; j++)
	{
		float temp = 0;
		float* pvalue = psf.ptr<float>(j);
		for (i = 0; i < pcol; i++)
		{
			temp += pvalue[i];
		}
		psfRowProj.at<float>(j) = temp;
	}
	CvSize tsize ;
	tsize.width = icol - 1;
	tsize.height = 1;

	copyMakeBorder(psfRowProj, psfRowProj, 0, tsize.height - psfRowProj.rows, 0, tsize.width - psfRowProj.cols, BORDER_CONSTANT, Scalar::all(0));
	fft2(psfRowProj, psfRowProj);

	Mat planes2[] = { Mat::zeros(tsize, CV_32F), Mat::zeros(tsize, CV_32F) };
	Mat temp1;
	split(psfRowProj, planes2);
	magnitude(planes2[0], planes2[1], planes2[0]);
	planes2[0] = planes2[0].mul(planes2[0]);
	dft(planes2[0], temp1, DFT_INVERSE + DFT_SCALE);
	split(psfRowProj, planes2);
	Mat beta[] = { Mat::zeros(1, icol, CV_32F), Mat::zeros(1, irow, CV_32F) };
	Mat temp2(beta[0], Rect(0, 0, icol - 1, 1));
	Mat temp3(beta[0], Rect(0, 0, 1, 1));
	Mat temp4(beta[0], Rect(icol - 1, 0, 1, 1));

	planes2[0].copyTo(temp2);
	temp3.copyTo(temp4);

	for (i = 0; i < pcol; i++)
	{
		float temp = 0;
		for (j = 0; j < prow; j++)
		{
			temp += psf.at<float>(j, i);
		}
		psfColProj.at<float>(i) = temp;
	}
	tsize.width = irow - 1;
	tsize.height = 1;

	copyMakeBorder(psfColProj, psfColProj, 0, tsize.height - psfColProj.rows, 0, tsize.width - psfColProj.cols, BORDER_CONSTANT, Scalar::all(0));
	fft2(psfColProj, psfColProj);

	Mat planes3[] = { Mat::zeros(tsize, CV_32F), Mat::zeros(tsize, CV_32F) };
	split(psfColProj, planes3);
	magnitude(planes2[0], planes2[1], planes2[0]);
	planes2[0] = planes2[0].mul(planes2[0]);
	dft(planes2[0], temp1, DFT_INVERSE + DFT_SCALE);
	dft(temp1, temp1, DFT_INVERSE + DFT_SCALE);
	split(psfColProj, planes3);
	Mat temp5(beta[1], Rect(0, 0, icol - 1, 1));
	Mat temp6(beta[1], Rect(0, 0, 1, 1));
	Mat temp7(beta[1], Rect(icol - 1, 0, 1, 1));

	planes3[0].copyTo(temp5);
	temp6.copyTo(temp7);

	Mat alpha(image.size(), CV_32F);
	beta[0] = 1 - beta[0];
	beta[1] = 1 - beta[1];
	transpose(beta[0], beta[0]);
	alpha = beta[0] * beta[1];
	Mat image1 = Mat_<float>(image);

	Mat J = alpha.mul(image1) + (1 - alpha).mul(blur);

	J = Mat_<uchar>(J);
	return J;
}