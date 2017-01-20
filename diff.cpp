#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <vector>

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

using namespace std;
using namespace cv;

int main(int argc, char** argv)
{
	Mat img1 = imread(argv[1],CV_LOAD_IMAGE_GRAYSCALE );
	Mat img2 = imread(argv[2],CV_LOAD_IMAGE_GRAYSCALE );

	Mat result(img1.size(),img1.type());

	absdiff(img1,img2,result);
	bitwise_not(result,result);

	imwrite("diff.pgm",result);
}