#include "basicoperation.h"
#include "estimatepar.h"
#include "restore.h"
#include "assessment.h"

void main()
{
	const string SRC_WIN_NAME = "Ô­Í¼";
	const string MED_WIN_NAME = "Ä£ºýÍ¼";
	const string DST_WIN_NAME = "¸´Ô­Í¼";
	const string SRC_NAME = "../../ImageRestore/images/flower.jpg";/*be_grdient.jpg*/


	Mat srcImage = imread(SRC_NAME, CV_LOAD_IMAGE_GRAYSCALE);
	Mat medImage, dstImage1, dstImage2, dstImage3, psf;
	imageInit(srcImage, srcImage);
	imshow("Ô­Í¼",srcImage);
	cout << "depth=" << srcImage.depth() << endl;
	fftShow(srcImage);
	//Mat A = (Mat_<float>(3, 3) << 3, 2, 1, 4.5, 5, 2, 4, 2, 3);
	//float d=A.at<float>(1,1);
	//float e = A.at<float>(1,2);
	//calWeight(srcImage, medImage);
	//cout << medImage<< endl;
	//normalize(medImage, medImage, 1, 0, NORM_MINMAX);
	//medImage=Mat_<uchar>(medImage);
	//genaratePsf(psf,15, 40);
//

//	filter2D(srcImage, medImage, -1, psf);
	//cepstral(medImage,dstImage1);
	//double  angle = directionDiffAngle(medImage);	
	//cout << " angle=" << angle << endl;
	//imshow("µ¹Æ×", dstImage1);
	//wnrFilter2(medImage, psf, dstImage1);

	//inverseFilter(medImage, psf, dstImage1);
//	deconvRL3(medImage, psf, dstImage1, 20);
	



	//double t1 = (double)getTickCount();
	//deconvRL4(medImage, psf, dstImage3, 20);
	//t1 = ((double)getTickCount() - t1) / getTickFrequency();

	//double t2 = (double)getTickCount();
	//deconvRL5(medImage, psf, dstImage2,20);
	//t2 = ((double)getTickCount() - t2) / getTickFrequency();
	////imwrite("lenal25a40.bmp", medImage);
	////imwrite("lenal15a40_wnr2.bmp", dstImage2);
	//deconvRL3(medImage, psf, dstImage3, 20);
	////double PSNR1 = calPSNR(dstImage2, srcImage);
	//double GMG1 = calGMG(srcImage);
	//double GMG2 = calGMG(medImage);
	//double GMG3 = calGMG(dstImage1);
	//double GMG4 = calGMG(dstImage2);
	//double con1 = calContrast(srcImage);
	//double con2 = calContrast(medImage);
	//double con3 = calContrast(dstImage1);
	//double con4 = calContrast(dstImage2);
	//double nsr1 = estimatek(srcImage);
	//double nsr2 = estimatek(medImage);
	//double nsr3 = estimatek(dstImage1);
	//double nsr4 = estimatek(dstImage2);


	////double PSNR3 = calPSNR(dstImage3, srcImage);

	////



	//
	////cout << "t=" <<t2 << endl;
	////
	//cout << "GMG1=" << GMG1 << endl;
	//cout << "GMG2=" << GMG2 << endl;
	//cout << "GMG3=" << GMG3 << endl;
	//cout << "GMG4=" << GMG4 << endl;
	//cout << "con1=" << con1 << endl;
	//cout << "con2=" << con2 << endl;
	//cout << "con3=" << con3 << endl;
	//cout << "con4=" << con4 << endl;
	//cout << "nsr1=" << nsr1 << endl;
	//cout << "nsr2=" << nsr2 << endl;
	//cout << "nsr3=" << nsr3 << endl;
	//cout << "nsr4=" << nsr4 << endl;

	////imshow("Ä£ºýºË", psf);
	//
	////namedWindow("Ä£ºýºË1", 0);

	//imshow(SRC_WIN_NAME, srcImage);
	//imshow(MED_WIN_NAME, medImage);
	//imshow("wnr+edgetaper", dstImage1);
	//imshow("RL5+edgetaper", dstImage2);
	//imshow("RL3+edgetaper", dstImage3);
//	imshow("extend", dstImage3);


	waitKey(0);
	destroyAllWindows();
}