#include"basicoperation.h"
#include"restore.h"
#include"assessment.h"
#include"estimatepar.h"

int main(int argc,char** argv)
{
	Mat srcImage;
	srcImage=imread(argv[1],0);
	imageInit(srcImage,srcImage);
	namedWindow("Display Image",CV_WINDOW_AUTOSIZE);
	imshow("Display Image",srcImage);

	fft2(srcImage,srcImage);
	fftshift(srcImage,srcImage);
	fftShow(srcImage);

	waitKey(0);

	return(0);
}
