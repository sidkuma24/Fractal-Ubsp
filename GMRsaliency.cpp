#include "GMRsaliency.h"

using namespace std;
typedef unsigned int UINT;

int isColor1 = 0;

Rect ROI(Mat &src);
void BbBShift(Mat& src, Mat& roiImage, Mat& dst);
void maxshift(Mat& src, Mat& roiImage, Mat& dst);
void bitslicing(Mat& src);

GMRsaliency::GMRsaliency()
{
	//
	spcount=200;
	compactness=20.0;
	alpha=0.99f;
	delta=0.1f;
	//
	spcounta=0;
}

GMRsaliency::~GMRsaliency()
{
}

Mat GMRsaliency::GetSup(const Mat &image)
{
	int width=image.cols;
	int height=image.rows;
	int sz = width*height;
    UINT *img=new UINT[sz*3];
	for(int c=0;c<3;c++)
	{
		for(int i=0;i<width;i++)
		{
			for(int j=0;j<height;j++)

				img[c*(width*height)+i*height+j]=saturate_cast<unsigned int>(image.at<Vec3b>(j,i)[2-c]);
		}
	}
	int* labels = new int[sz];

	SLIC slic;
	slic.DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(img, height, width, labels, spcounta, spcount, compactness);	


	Mat supLab(image.size(),CV_16U);
	for(int i=0;i<supLab.rows;i++)
	{
		for(int j=0;j<supLab.cols;j++)
		{
			supLab.at<ushort>(i,j)=labels[i+j*supLab.rows];

		}
	}


	if(labels) delete [] labels;
	if(img) delete []img;

	return supLab;

}

Mat GMRsaliency::GetAdjLoop(const Mat &supLab)
{
	Mat adj(Size(spcounta,spcounta),CV_16U,Scalar(0));
	for(int i=0;i<supLab.rows-1;i++)
	{
		for(int j=0;j<supLab.cols-1;j++)
		{
			if(supLab.at<ushort>(i,j)!=supLab.at<ushort>(i+1,j))
			{
				adj.at<ushort>(supLab.at<ushort>(i,j),supLab.at<ushort>(i+1,j))=1;
				adj.at<ushort>(supLab.at<ushort>(i+1,j),supLab.at<ushort>(i,j))=1;
			}
			if(supLab.at<ushort>(i,j)!=supLab.at<ushort>(i,j+1))
			{
				adj.at<ushort>(supLab.at<ushort>(i,j),supLab.at<ushort>(i,j+1))=1;
				adj.at<ushort>(supLab.at<ushort>(i,j+1),supLab.at<ushort>(i,j))=1;
			}
			if(supLab.at<ushort>(i,j)!=supLab.at<ushort>(i+1,j+1))
			{
				adj.at<ushort>(supLab.at<ushort>(i,j),supLab.at<ushort>(i+1,j+1))=1;
				adj.at<ushort>(supLab.at<ushort>(i+1,j+1),supLab.at<ushort>(i,j))=1;
			}
			if(supLab.at<ushort>(i+1,j)!=supLab.at<ushort>(i,j+1))
			{
				adj.at<ushort>(supLab.at<ushort>(i+1,j),supLab.at<ushort>(i,j+1))=1;
				adj.at<ushort>(supLab.at<ushort>(i,j+1),supLab.at<ushort>(i+1,j))=1;
			}
		}
	}
	vector<ushort>bd;
	vector<ushort>::iterator result;
	for(int i=0;i<supLab.cols;i++)
	{
		result=find(bd.begin(),bd.end(),supLab.at<ushort>(0,i));
		if(result==bd.end())
			bd.push_back(supLab.at<ushort>(0,i));
	}
	for(int i=0;i<supLab.cols;i++)
	{
		result=find(bd.begin(),bd.end(),supLab.at<ushort>(supLab.rows-1,i));
		if(result==bd.end())
			bd.push_back(supLab.at<ushort>(supLab.rows-1,i));
	}
	for(int i=0;i<supLab.rows;i++)
	{
		result=find(bd.begin(),bd.end(),supLab.at<ushort>(i,0));
		if(result==bd.end())
			bd.push_back(supLab.at<ushort>(i,0));
	}
	for(int i=0;i<supLab.rows;i++)
	{
		result=find(bd.begin(),bd.end(),supLab.at<ushort>(i,supLab.cols-1));
		if(result==bd.end())
			bd.push_back(supLab.at<ushort>(i,supLab.cols-1));
	}
	vector<ushort>::iterator bdi=bd.begin();
	vector<ushort>::iterator bdj;
	for(;bdi!=bd.end();bdi++)
	{
		for(bdj=bdi+1;bdj!=bd.end();bdj++)
		{
			adj.at<ushort>(*bdi,*bdj)=1;
			adj.at<ushort>(*bdj,*bdi)=1;
		}
	}

	return adj;

}

Mat GMRsaliency::GetWeight(const Mat &img,const Mat &supLab,const Mat &adj)
{
	vector<float> supL(spcounta,0);
	vector<float> supa(spcounta,0);
	vector<float> supb(spcounta,0);
	vector<float> pcount(spcounta,0);
	for(int i=0;i<img.rows;i++)
	{
		const float *labp=img.ptr<float>(i);
		for(int j=0;j<img.cols;j++,labp+=3)
		{
			supL[supLab.at<ushort>(i,j)]+=labp[0];
			supa[supLab.at<ushort>(i,j)]+=labp[1];
			supb[supLab.at<ushort>(i,j)]+=labp[2];
			pcount[supLab.at<ushort>(i,j)]+=1.;
		}
	}
	for(int i=0;i<spcounta;i++)
	{
		supL[i]/=pcount[i];
		supa[i]/=pcount[i];
		supb[i]/=pcount[i];
	}
	Mat w(adj.size(),CV_32F,Scalar(-1));
	float minw=(float)numeric_limits<float>::max(),maxw=(float)numeric_limits<float>::min();
	for(int i=0;i<spcounta;i++)
	{
		for(int j=0;j<spcounta;j++)
		{
			if(adj.at<ushort>(i,j)==1)
			{
				float dist=sqrt(pow((supL[i]-supL[j]),2)+pow((supa[i]-supa[j]),2)+pow((supb[i]-supb[j]),2));
				w.at<float>(i,j)=dist;
				if(minw>dist)
					minw=dist;
				if(maxw<dist)
					maxw=dist;
				for(int k=0;k<spcounta;k++)
				{
					if(adj.at<ushort>(j,k)==1&&k!=i)
					{
						float dist=sqrt(pow((supL[i]-supL[k]),2)+pow((supa[i]-supa[k]),2)+pow((supb[i]-supb[k]),2));
						w.at<float>(i,k)=dist;
						if(minw>dist)
							minw=dist;
						if(maxw<dist)
							maxw=dist;
					}
						
				}
			}
		}
	}
	for(int i=0;i<spcounta;i++)
	{
		for(int j=0;j<spcounta;j++)
		{
			if(w.at<float>(i,j)>-1)
				w.at<float>(i,j)=exp(-(w.at<float>(i,j)-minw)/((maxw-minw)*delta));
			else
				w.at<float>(i,j)=0;
		}
	}
	return w;
}

Mat GMRsaliency::GetOptAff(const Mat &W)
{

	Mat dd(Size(W.rows,1),CV_32F);
	reduce(W,dd,1,CV_REDUCE_SUM);
	Mat D(W.size(),CV_32F);
	D=Mat::diag(dd);
	Mat optAff(W.size(),CV_32F);
	optAff=(D-alpha*W);
	optAff=optAff.inv();
	Mat B=Mat::ones(optAff.size(),CV_32F)-Mat::eye(optAff.size(),CV_32F);
	optAff=optAff.mul(B);
	return optAff;

}

Mat GMRsaliency::GetBdQuery(const Mat &supLab,int type)
{
	Mat y(Size(1,spcounta),CV_32F,Scalar(0));
	switch(type)
	{
	case 1:
		for(int i=0;i<supLab.cols;i++)
			y.at<float>(supLab.at<ushort>(0,i))=1;
		break;
	case 2:
		for(int i=0;i<supLab.cols;i++)
			y.at<float>(supLab.at<ushort>(supLab.rows-1,i))=1;
		break;
	case 3:
		for(int i=0;i<supLab.rows;i++)
			y.at<float>(supLab.at<ushort>(i,0))=1;
		break;
	case 4:
		for(int i=0;i<supLab.rows;i++)
			y.at<float>(supLab.at<ushort>(i,supLab.cols-1))=1;
		break;
	default:
		printf("error");
	}
	return y;
		
}

Mat GMRsaliency::RemoveFrame(const Mat &img,int *wcut)
{
	double thr=0.6;
	Mat grayimg;
	cvtColor(img,grayimg,CV_BGR2GRAY);
	Mat edge;
	Canny(grayimg,edge,150*0.4,150);

	int flagt=0;
	int flagd=0;
	int flagr=0;
	int flagl=0;
	int t=0;
	int d=0;
	int l=0;
	int r=0;

	int m=grayimg.rows;
	int n=grayimg.cols;

	int i=0;
	while(i<30)
	{
		float pbt=(float)mean(edge(Range(i,i+1),Range::all()))[0];
		float pbd=(float)mean(edge(Range(m-i-1,m-i),Range::all()))[0];
		float pbl=(float)mean(edge(Range::all(),Range(i,i+1)))[0];
		float pbr=(float)mean(edge(Range::all(),Range(n-i-1,n-i)))[0];
		if(pbt/255>thr)
		{
			t=i;
			flagt=1;
		}
		if(pbd/255>thr)
		{
			d=i;
			flagd=1;
		}
		if(pbl/255>thr)
		{
			l=i;
			flagl=1;
		}
		if(pbr/255>thr)
		{
			r=i;
			flagr=1;
		}
		i++;
	}
	int flagrm=flagt+flagd+flagl+flagr;
	Mat outimg;
	if(flagrm>1)
	{
		int maxwidth;
		maxwidth=(t>d)?t:d;
		maxwidth=(maxwidth>l)?maxwidth:l;
		maxwidth=(maxwidth>r)?maxwidth:r;
		if(t==0)
			t=maxwidth;
		if(d==0)
			d=maxwidth;
		if(l==0)
			l=maxwidth;
		if(r==0)
			r=maxwidth;
		outimg=img(Range(t,m-d),Range(l,n-r));
		wcut[0]=m;
		wcut[1]=n;
		wcut[2]=t;
		wcut[3]=m-d;
		wcut[4]=l;
		wcut[5]=n-r;
	}
	else
	{
		wcut[0]=m;
		wcut[1]=n;
		wcut[2]=0;
		wcut[3]=m;
		wcut[4]=0;
		wcut[5]=n;
		outimg=img;
	}
	return outimg;

}

void quicksort(float *x,int first,int last){
    float temp;
    int j, i, pivot;
    

     if(first<last){
         pivot=first;
         i=first;
         j=last;

         while(i<j){
             while(x[i]<=x[pivot]&&i<last)
                 i++;
             while(x[j]>x[pivot])
                 j--;
             if(i<j){
                 temp=x[i];
                  x[i]=x[j];
                  x[j]=temp;
             }
         }

         temp=x[pivot];
         x[pivot]=x[j];
         x[j]=temp;
         quicksort(x,first,j-1);
         quicksort(x,j+1,last);

    }
}

 float median(float * array, int size)
 {
 	quicksort(array,0,size-1);

 	int middle;
 	float median;	
    middle = size / 2.0;

    if(size % 2 == 0){
    	median = (array[middle] + array[middle + 1] )/2.0;
    }
    else{
    	median = array[middle];
    }

    return median;

 }

Mat GMRsaliency::GetSal(Mat &img)
{

	int wcut[6];

	img=RemoveFrame(img,wcut);
	//cout << "img" << img << endl;
     //imshow("image",img); waitKey(0);
	Mat suplabel(img.size(),CV_16U);
	suplabel=GetSup(img);
   
    Mat adj(Size(spcounta,spcounta),CV_16U);
	adj=GetAdjLoop(suplabel);
   // imshow("image",adj); waitKey(0);

	Mat tImg;
	img.convertTo(tImg,CV_32FC3,1.0/255);
	cvtColor(tImg, tImg, CV_BGR2Lab);
    
	Mat W(adj.size(),CV_32F);
	W=GetWeight(tImg,suplabel,adj);
   
	Mat optAff(W.size(),CV_32F);
	optAff=GetOptAff(W);
   // cout << "optAff" << optAff << endl;
	Mat salt,sald,sall,salr;
	Mat bdy;

	bdy=GetBdQuery(suplabel,1);
	salt=optAff*bdy;
	normalize(salt, salt, 0, 1, NORM_MINMAX);
	salt=1-salt;

	bdy=GetBdQuery(suplabel,2);
	sald=optAff*bdy;
	normalize(sald, sald, 0, 1, NORM_MINMAX);
	sald=1-sald;

	bdy=GetBdQuery(suplabel,3);
	sall=optAff*bdy;
	normalize(sall, sall, 0, 1, NORM_MINMAX);
	sall=1-sall;

	bdy=GetBdQuery(suplabel,4);
	salr=optAff*bdy;
	normalize(salr, salr, 0, 1, NORM_MINMAX);
	salr=1-salr;

	Mat salb;
	salb=salt;
	salb=salb.mul(sald);
	salb=salb.mul(sall);
	salb=salb.mul(salr);

	double thr=mean(salb)[0];
	Mat fgy;
	threshold(salb,fgy,thr,1,THRESH_BINARY);
    //imshow("fgy",fgy); waitKey(0);
    Mat salf;
	salf=optAff*fgy;
	
	//cout << "salf rows" <<salf.rows << endl;
	//cout << "salf cols" <<salf.cols << endl;
    cout << salf << endl;
	float salf_array[salf.rows];
	for(int r=0; r < salf.rows; r++){
		salf_array[r] = salf.at<float>(r,1);
		//cout << "array:" << salf_array[r] << endl;
	}
    float med = median(salf_array,salf.rows);
    //cout << "Median :" << med << endl; 
    for(int r=0; r < salf.rows; r++){
		if(salf.at<float>(r,0) < med){
			salf.at<float>(r,0) = 0;
		}else{
			salf.at<float>(r,0) = 1;
		}
		//cout << "array:" << salf_array[r] << endl;
	}


 
	Mat salMap(img.size(),CV_32F);
	for(int i=0;i<salMap.rows;i++)
	{
		for(int j=0;j<salMap.cols;j++)
		{
			salMap.at<float>(i,j)=salf.at<float>(suplabel.at<ushort>(i,j));
		}
	}
    
	normalize(salMap, salMap, 0, 1, NORM_MINMAX);
//	imshow("sal map",salMap); waitKey(0);
  //  cout << "sal final "<< salf << endl ;
	Mat outMap(Size(wcut[1],wcut[0]),CV_32F,Scalar(0));
	Mat subMap=outMap(Range(wcut[2],wcut[3]),Range(wcut[4],wcut[5]));
	salMap.convertTo(subMap,subMap.type());
	Mat dst;
	if(isColor1){
		dst = (img.size(),CV_16UC3);
	}else{
		dst = (img.size(),CV_16U);
	}
	maxshift(img,outMap,dst);
	Mat output(img.size(),CV_8U);
	dst.convertTo(output,CV_8U);
	imwrite("sal.png",subMap);
//	BbBShift(img, outMap, dst);
//	bitslicing(dst);
	return outMap;

}


void maxshift(Mat& src, Mat& roiImage, Mat& dst)
{
	cvtColor(src,src,CV_RGB2GRAY);
	int shift = 4;
	Mat shifted;
	if(isColor1){
		shifted = Mat(src.size(),CV_16UC3);
	}else{
		shifted = Mat(src.size(),CV_16U);
		for(int r =0; r < shifted.rows; ++r){
			for(int c=0; c < shifted.cols; ++c){
				shifted.at<unsigned short>(r,c) = src.at<uchar>(r,c);
				if(roiImage.at<uchar>(r,c) == 0){
					shifted.at<unsigned short>(r,c) = shifted.at<unsigned short>(r,c) >> shift; 
				}
			}
		}
	}
	shifted.copyTo(dst);
}


void BbBShift(Mat& src, Mat& roiImage, Mat& dst)
{
	Mat gray(src.size(),CV_8U);
	cvtColor(src,gray,CV_RGB2GRAY);
	roiImage.convertTo(roiImage,CV_8U);

	int value =0;
	unsigned short bp_mask = 61680;
	
    Mat src_shifted(src.size(),CV_16U);
    if(isColor1){
    	src_shifted.convertTo(src_shifted,CV_8UC3);
    }
    Mat ROI, nonROI;
    if(isColor1){
	}else{
      	ROI = Mat(src.size(),CV_8U);
      	nonROI = Mat(src.size(),CV_8U);
      	for(int r=0; r < src.rows; r++){
      		for(int c=0; c < src.cols; c++){
      			if(roiImage.at<uchar>(r,c) == 1 ){
      				ROI.at<uchar>(r,c) = gray.at<uchar>(r,c);
       			}else{
       				nonROI.at<uchar>(r,c) = gray.at<uchar>(r,c);
       			}
      		}
    	}
    }
    imshow("gray",gray); waitKey(0);
    imshow("Non ROI ",nonROI); waitKey(0);
      unsigned short num = 32768;
      unsigned char bitPos1 = 128;
      unsigned char bitPos2 = 128; 
      int roi_count  = 0;
      int non_roi_count = 0;
      for(int i=0; i < 16; i++){
      	if((bp_mask & num)){
      		roi_count++;
      		cout << "ROI count =" << roi_count << endl;
      		for(int r=0; r < src.rows; r++){
      			for(int c=0; c < src.cols; c++){
      		
		      		src_shifted.at<unsigned short>(r,c) = src_shifted.at<unsigned short>(r,c) | (ROI.at<uchar>(r,c) & bitPos1);
		      		src_shifted.at<unsigned short>(r,c) << 8 - non_roi_count ;
		      	}
		    }

		    bitPos1 = bitPos1 >> 1;

		}else{
			cout << "Non roi count = " << non_roi_count << endl;
			non_roi_count++;
			for(int r=0; r < src.rows; r++){
      			for(int c=0; c < src.cols; c++){
      				if(roiImage.at<uchar>(r,c) == 0){
		      		src_shifted.at<unsigned short>(r,c) = src_shifted.at<unsigned short>(r,c) | (nonROI.at<uchar>(r,c) & bitPos2);
		      		src_shifted.at<unsigned short>(r,c) << 8 - roi_count;

		      		src_shifted.at<unsigned short>(r,c) = src_shifted.at<unsigned short>(r,c) & 65472;
		      	}
		      	}
		    }
		    bitPos2 = bitPos2 >> 1;
		}
		bp_mask = bp_mask << 1;

	}
	cout << "ROI count =" << roi_count << endl;
	cout << "Non roi count = " << non_roi_count << endl;
	// int num = 128;
	src_shifted.copyTo(dst);

	// cout <<"128 after 6 right shifts= " << (128 >> 6) << endl;
	//imwrite("final.png",final);
}


void bitslicing(Mat& src)
{
	Mat gray(src.size(),src.type());
	if(src.channels() != 1)
		cvtColor(src,gray,CV_RGB2GRAY);
	else
		src.copyTo(gray);
    gray.convertTo(gray,CV_8U);
    imshow("out",gray);waitKey(0);
	Mat out(src.size(),CV_8U);
	unsigned int num = 1;
    int lowest_bit = 8;
    
    for(int r=0 ; r < src.rows; ++r){
    	for(int c = 0;  c < src.cols; ++c){
    		out.at<uchar>(r,c) = gray.at<unsigned short>(r,c) >> 8;
		}
	}

	imwrite("outImage.png",gray);

	
}
