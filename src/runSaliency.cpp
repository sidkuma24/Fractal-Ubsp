/*
demo for paper "Saliency Detection via Graph-based Manifold Ranking"
by Chuan Yang, Lihe Zhang, Huchuan Lu, Xiang Ruan, and Ming-Hsuan Yang
in CVPR13.
written by Chuan Yang
email: ycscience86@gmail.com
date: 2013.5.7
*/

#include <cstring>
#include <stdio.h>
#include <time.h>
#include "GMRsaliency.h"

void GetImageNames(char* filepath,vector<string> &imnames);

int main(int argc,char *argv[])
{
	// char filepath[MAX_PATH+1];
	// sprintf(filepath,"%s*%s",argv[1],argv[2]);
	// vector<string> imnames;
	// GetImageNames(filepath,imnames);
	Mat sal,img;
	//vector<string>::iterator itf=imnames.begin();
	int count=1;
	clock_t nTimeStart;      
    clock_t nTimeStop;       
    nTimeStart = clock();

	img=imread(argv[1]);
	//bit_mask = argv[2];
	//bit_position = atoi(argv[2]);
	// if(img.channels() == 3)
	// 	isColor = 1;
	GMRsaliency GMRsal;
	sal=GMRsal.GetSal(img);
	
	//sprintf(salname,"%s_our.png",imname);
	imwrite("sal_image.png",sal*255);

	nTimeStop = clock();
	cout <<"the average running time:"<<(double)(nTimeStop - nTimeStart)/CLOCKS_PER_SEC/(count-1)<<"seconds"<< endl;

	system("pause");
	return 0;
}
// void GetImageNames(char* filepath,vector<string> &imnames)
// {
// 	imnames.clear();
// 	WIN32_FIND_DATA f;
// 	HANDLE h = FindFirstFile(filepath , &f);
// 	if(h != INVALID_HANDLE_VALUE)
// 	{
// 		do
// 		{
// 			imnames.push_back(f.cFileName);
// 		} while(FindNextFile(h, &f));

// 	}
// 	FindClose(h);	
// }
