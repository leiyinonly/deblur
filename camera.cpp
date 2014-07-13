#include <cv.h>
#include <opencv/highgui.h>
#include<stdio.h>
#include<opencv2/opencv.hpp>
using namespace cv;

int main(int argc, char *argv[])
{
      CvCapture* pCapture = cvCreateCameraCapture(0);
      cvNamedWindow("Video", 1);

      while(1)
      {
          IplImage* pFrame=cvQueryFrame( pCapture );
          if(!pFrame)break;
          cvShowImage("Video",pFrame);
          char c=cvWaitKey(33);
          if(c==27)break;
      }
      cvReleaseCapture(&pCapture);
      cvDestroyWindow("Video");
      return 0;
}
