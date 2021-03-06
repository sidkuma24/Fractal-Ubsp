/******************************************************************************
 ==============================================================================

             '`
            '--`        Mars -- A quadtree based fractal image
           '`  '`       coder/decoder. Version 1.0 (10/28/1998) 
          '--`'--`      Copyright (C) 1998 Mario Polvere 
         '`      '`     University of Salerno Italy 
        '--`    '--`    
       '`  '`  '`  '`
      '--`'--`'--`'--`

      This program is free software; you can redistribute it and/or modify
      it under the terms of the GNU General Public License as published by
      the Free Software Foundation; either version 2 of the License, or
      (at your option) any later version.

      This program is distributed in the hope that it will be useful,
      but WITHOUT ANY WARRANTY; without even the implied warranty of
      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
      GNU General Public License for more details.

      You should have received a copy of the GNU General Public License
      along with this program; if not, write to the Free Software
      Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 
      02111-1307, USA

      The author can be contacted via e-mail:
                               marpol@iname.com
                               teimars@erifs1.ericsson.se

 ==============================================================================
 ******************************************************************************/



#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <vector>

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include "def.h"
#define EXTERN
#include "globals.h"
#include "triangle.h"
#include "prot.h"

 using namespace cv;
 
int main(int argc, char **argv)
{
  int i,j,k,max;
  char *p,*pp;
  int int_max_alfa;
  clock_t t;
  double time_taken;


 
  printf("%s",triangle_3);

  getopt_enc(argc,argv);

  // pp = "";
  // p = filein;
  // while(*p) {        Find the last dot in the file name 
  //   if(*p == '.') 
  //      pp = p + 1;
  //   p++;
  // }

  // if(!strcmp(pp,"pgm"))
  //    readimage_pgm(filein,&image_width,&image_height);  
  // else
  //    readimage_raw(filein);  

  Mat input_image = imread(filein);
  int w = input_image.cols;
  int h = input_image.rows;   

  image_width = w;
  image_height = h;

  std::vector<Mat> yuvChannels(3);
  if(input_image.channels() == 3){
    isColor = 1;
    
    cvtColor(input_image, input_image, CV_BGR2YUV);
    split(input_image, yuvChannels);
    matrix_allocate(image, w+1, h+1, PIXEL);
    matrix_allocate(image_uch, w+1, h+1, PIXEL);
    matrix_allocate(image_vch, w+1, h+1, PIXEL);
    for (int iii = 0; iii < w; iii++) {
      for (int jjj = 0; jjj < h; jjj++) {
        image[jjj][iii] = yuvChannels[0].at<uchar>(jjj, iii);
        image_uch[jjj][iii] = yuvChannels[1].at<uchar>(jjj, iii);
        image_vch[jjj][iii] = yuvChannels[2].at<uchar>(jjj, iii);
      }
    }
  }
  else{
    cvtColor(input_image,input_image,CV_RGB2GRAY);
    matrix_allocate(image, w, h, PIXEL);
    for (int iii = 0; iii < w; iii++) {
      for (int jjj = 0; jjj < h; jjj++) {
        image[jjj][iii] = input_image.at<uchar>(jjj, iii);
      }
    }
   }


  max = image_height;
  if(image_width > image_height ) 
    max = image_width; 

  virtual_size = 1 << (int) ceil(log((double) max) / log(2.0));

  matrix_allocate(contract,1+image_width / 2,1+image_height / 2,double)
  matrix_allocate(qtt,virtual_size,virtual_size,PIXEL)
  matrix_allocate(hv,virtual_size+1,virtual_size+1,PIXEL)
  matrix_allocate(range_tmp,64,64,double)
  matrix_allocate(flip_range,64,64,double)
  matrix_allocate(range,64,64,double)
  t = clock();
  contraction(contract,image,0,0);
  


  switch(method) {
    case MassCenter  :  
      ComputeAverageFactorMc();
      for(k=0;k<8;k++) {
         matrix_allocate(class_polar[k],n_p_class,n_p_class,struct c *)
      }

      for(k=0;k<8;k++)
      for(h=0;h<n_p_class;h++)
      for(i=0;i<n_p_class;i++)
         class_polar[k][h][i] = NULL;

      Indexing = MassCenterIndexing; 
      Coding = MassCenterCoding;
      printf(" Speed-up method: MassCenter\n\n");
      break;
    case SaupeFisher :  
      ComputeFeatVectDimSaupe();
      Indexing = Saupe_FisherIndexing;
      Coding = Saupe_FisherCoding;
      printf(" Speed-up method: Saupe-Fisher\n\n");
      break;
    case Saupe :  
      ComputeFeatVectDimSaupe();
      Indexing = SaupeIndexing;
      Coding = SaupeCoding;
      printf(" Speed-up method: Saupe\n\n");
      break;
     case CoVar :
      // for(i=(int) rint(log((double) 2) / log(2.0)); 
      //                i<= (int) rint(log((double) max_size) / log(2.0)); i++) {
      //  findMaxStd((int) rint(pow(2.0,(double)i)),i);
      // }  
      // for(int i=0; i < 8; ++i){
      //   if(final_max_std < max_std_arr[i]) final_max_std = max_std_arr[i];
      // }
      // final_max_std = (int)floor(final_max_std);
      final_max_std = 255;
      // printf("Max STD = %f\n",final_max_std);
      // for(int i=0;i<8;++i){
      //   for(int j=0; j<3;++j){
      //     class_std[i][j] = (struct c **) malloc(sizeof(struct c*) * (final_max_std+5));
      //   }
      // }
      for(k=0;k<8;k++)
      for(h=0;h<3;h++)
      for(i=0;i<final_max_std;++i)
          class_std[k][h][i] = NULL;

      // for(i=0;i<24;i++)
          // class_basic[k][h][i] = NULL;

      Indexing = COVIndexing;
      Coding =  COVCoding;
      printf(" Speed-up method: Coeff of Variance  \n\n");
      break;

      case CoVar1 :
      // for(i=(int) rint(log((double) 2) / log(2.0)); 
      //                i<= (int) rint(log((double) max_size) / log(2.0)); i++) {
      //  findMaxStd((int) rint(pow(2.0,(double)i)),i);
      // }  
      // for(int i=0; i < 8; ++i){
      //   if(final_max_std < max_std_arr[i]) final_max_std = max_std_arr[i];
      // }
      // final_max_std = (int)floor(final_max_std);
      final_max_std = 255;
      // printf("Max STD = %f\n",final_max_std);
     
      for(k=0;k<8;k++)
      for(h=0;h<3;h++)
      for(i=0;i<final_max_std;++i)
          class_std[k][h][i] = NULL;

      // for(i=0;i<24;i++)
          // class_basic[k][h][i] = NULL;

      Indexing = COVIndexing;
      Coding =  CovClass_AdaptiveSearch_FisherCoding;
      // Coding =  COVCoding;
      printf(" Speed-up method: Coeff of Variance with adaptive search  \n\n");
      break;

       case CoVar2 :
      // for(i=(int) rint(log((double) 2) / log(2.0)); 
      //                i<= (int) rint(log((double) max_size) / log(2.0)); i++) {
      //  findMaxStd((int) rint(pow(2.0,(double)i)),i);
      // }  
      // for(int i=0; i < 8; ++i){
      //   if(final_max_std < max_std_arr[i]) final_max_std = max_std_arr[i];
      // }
      // final_max_std = (int)floor(final_max_std);
      final_max_std = 255;
      // printf("Max STD = %f\n",final_max_std);
      // for(int i=0;i<8;++i){
      //   for(int j=0; j<3;++j){
      //     class_std[i][j] = (struct c **) malloc(sizeof(struct c*) * (final_max_std+5));
      //   }
      // }
      for(k=0;k<8;k++)
      for(h=0;h<3;h++)
      for(i=0;i<final_max_std;++i)
          class_std[k][h][i] = NULL;

      // for(i=0;i<24;i++)
          // class_basic[k][h][i] = NULL;

      Indexing = COVIndexing;
      Coding =  CovClass_AdaptiveSearch_FisherCoding;
      printf(" Speed-up method: Coeff of Variance with adaptive search and adaptive range based paritioning  \n\n");
      break;

    case Fisher_std :
      for(i=(int) rint(log((double) 2) / log(2.0)); 
                     i<= (int) rint(log((double) max_size) / log(2.0)); i++) {
       findMaxStd((int) rint(pow(2.0,(double)i)),i);
      }  
      for(int i=0; i < 8; ++i){
        if(final_max_std < max_std_arr[i]) final_max_std = max_std_arr[i];
      }
      final_max_std = 255;
      // printf("Max STD = %f\n",final_max_std);
      // for(int i=0;i<8;++i){
      //   for(int j=0; j<3;++j){
      //     class_std[i][j] = (struct c **) malloc(sizeof(struct c*) * (final_max_std+5));
      //   }
      // }
      for(k=0;k<8;k++)
      for(h=0;h<3;h++)
      for(i=0;i<final_max_std;++i)
          class_std[k][h][i] = NULL;

      // for(i=0;i<24;i++)
          // class_basic[k][h][i] = NULL;

      Indexing = STDIndexing;
      Coding =  STDCoding;
      printf(" Speed-up method: STD \n\n");
      break;
    case Entropy_based :
      for(i=(int) rint(log((double) 2) / log(2.0)); 
                     i<= (int) rint(log((double) max_size) / log(2.0)); i++) {
       findMaxEnt((int) (rint(pow(2.0,(double)i))/2),i);
      }  
      for(int i=0; i < 8; ++i){
        if(final_max_ent > max_ent_arr[i]) final_max_ent = max_ent_arr[i];
      }
      final_max_ent_q = (int)floor(final_max_ent);
      printf("Max ENT = %d\n",final_max_ent_q);
      for(int i=0;i<8;++i){
      for(int j=0; j<3;++j){
        class_entropy[i][j] = (struct c **) malloc(sizeof(struct c*) * final_max_ent_q);
        }
      }
      for(k=0;k<8;k++)
      for(h=0;h<3;h++)
      for(i=0;i<final_max_ent_q;++i)
          class_entropy[k][h][i] = NULL;

      // for(i=0;i<24;i++)
          // class_basic[k][h][i] = NULL;

      Indexing = EntropyIndexing;
    Coding =  EntropyCoding;
      printf(" Speed-up method: Entropy based classification\n\n");
      break;

    case Fisher :  
      for(k=0;k<8;k++)
      for(h=0;h<3;h++)
      for(i=0;i<24;i++)
          class_fisher[k][h][i] = NULL;

      Indexing = FisherIndexing;
      Coding = FisherCoding;
      printf(" Speed-up method: Fisher Coding\n\n");
      break;
      
    case BasicFIC :  
      for(k=0;k<8;k++)
      for(h=0;h<3;h++)
          class_basicFIC[k][h] = NULL;

      Indexing = BasicFIC_Indexing;
      Coding = BasicFIC_Coding;
      printf(" Speed-up method: Basic FIC\n\n");
      break;
    
    case Nonlinear_Fisher :  
      for(k=0;k<8;k++)
      for(h=0;h<3;h++)
      for(i=0;i<24;i++)
          class_fisher[k][h][i] = NULL;

      Indexing = FisherIndexing;
      Nonlinear_Coding = Nonlinear_FisherCoding;
      printf(" Speed-up method: Fisher nonlinear\n\n");
      break;
    case LumInv_Fisher :  
      for(k=0;k<8;k++)
      for(h=0;h<3;h++)
      for(i=0;i<24;i++)
          class_fisher[k][h][i] = NULL;

      Indexing = FisherIndexing;
      testing_Coding = testing_FisherCoding;
      printf(" Speed-up method:Luminance Invariant Fisher \n\n");
      break;
    case testing_Fisher :  
      for(k=0;k<8;k++)
      for(h=0;h<3;h++)
      for(i=0;i<24;i++)
          class_fisher[k][h][i] = NULL;

      Indexing = FisherIndexing;
      testing_Coding = testing_FisherCoding;
      printf(" Speed-up method:Luminance Invariant Fisher \n\n");
      break;
    case Hurtgen :  
      for(k=0;k<8;k++)
      for(h=0;h<16;h++)
      for(i=0;i<24;i++)
          class_hurtgen[k][h][i] = NULL;
      Indexing = HurtgenIndexing;
      Coding = HurtgenCoding;
      printf(" Speed-up method: Hurtgen\n\n");
      break;
    case McSaupe :  
      ComputeAverageFactorMc();
      ComputeFeatVectDimSaupe();

      for(k=0;k<8;k++) {
         matrix_allocate(class_polar[k],n_p_class,n_p_class,struct c *)
         matrix_allocate(class_polar_saupe[k],n_p_class,n_p_class,kdtree *)
         matrix_allocate(item_per_class[k],n_p_class,n_p_class, int )
         matrix_allocate(c_book[k],n_p_class,n_p_class, struct code_book *)
         matrix_allocate(f_vect[k],n_p_class,n_p_class, float ** )
      }

      for(k=0;k<8;k++)
      for(h=0;h<n_p_class;h++)
      for(i=0;i<n_p_class;i++) {
          class_polar[k][h][i] = NULL; 
          class_polar_saupe[k][h][i] = NULL; 
          item_per_class[k][h][i] = 0; 
      }
      Indexing = Mc_SaupeIndexing;
      Coding = Mc_SaupeCoding;
      printf(" Speed-up method: Mc-Saupe\n\n");
      break;
 
    /*
     * Fisher Classification with modified quadtree partition 
     */
    case Fisher_adaptiveQuadtree:
      int l;
      for(k=0;k<64;k++)
      for(l=0;l<64;l++)
      for(h=0;h<24;h++)
          adaptive_fisher_class[k][l][h] = NULL;

      adaptiveIndexing = adaptiveFisherIndexing_2;
      adaptiveCoding = adaptiveFisherCoding;
      printf(" Speed-up method: Fisher with Adaptive Quadtree\n\n");
    break;

      case modified_Fisher :  
      for(k=0;k<8;k++)
      for(h=0;h<3;h++)
      for(i=0;i<24;i++)
          class_fisher[k][h][i] = NULL;

      Indexing = FisherIndexing;
      Coding = modified_FisherCoding;
      printf(" Speed-up method: Modified Fisher Coding\n\n");
      break;

       case modified_Fisher1 :  
      for(k=0;k<8;k++)
      for(h=0;h<3;h++)
      for(i=0;i<24;i++)
          class_fisher[k][h][i] = NULL;

      Indexing = FisherIndexing;
      Coding = modified_FisherCoding_1;
      printf(" Speed-up method: Modified Fisher Coding with HV patition\n\n");
      break;

       case modified_Fisher2 :  
      for(k=0;k<8;k++)
      for(h=0;h<3;h++)
      for(i=0;i<24;i++)
          class_fisher[k][h][i] = NULL;

      Indexing =  FisherIndexing;
      Coding = modified_FisherCoding_2;
      printf(" Speed-up method: Modified Fisher Coding with HV patition and adaptive search\n\n");
      break;

       case modified_Fisher3 :  
      for(k=0;k<8;k++)
      for(h=0;h<3;h++)
      for(i=0;i<24;i++)
          class_fisher[k][h][i] = NULL;

      Indexing = FisherIndexing_domainSort;
      Coding = modified_FisherCoding_2;
      printf(" Speed-up method: Modified Fisher Coding with HV patition and adaptive search\n"
             " and domain sorting                                                        \n\n");
      break;
    /****************************************************************/
    case Fisher_HV:
      
      for(k=0;k<64;k++)
      for(l=0;l< 64 ;l++)
      for(h=0;h<24;h++){
          adaptive_fisher_class[k][l][h] = NULL;

      }
      HV_Indexing = HV_FisherIndexing;
      HV_Coding = HV_FisherCoding;
      printf(" Speed-up method: Fisher with HV parition\n\n");
    break;
   /*
    case Your_method :                               If you want to try a new 
      classification = YourMethodIndexing;           method you need just to
      coding = YourMethodCoding;                     write an indexing and a
      printf(" Speed-up method: YourMethod\n\n");    coding function and call
      break;                                         them here 
   */
  } 

  // contraction(contract,image,0,0);

  switch(partition_type){

    case Quadtree:
       for(i=(int) rint(log((double) 2) / log(2.0)); 
                     i<= (int) rint(log((double) max_size) / log(2.0)); i++) {
          Indexing((int) rint(pow(2.0,(double)i)),i);
        }
    break;

    case AdaptiveQuadtree:
        
        for(i= 1 << MIN_ADAP_D_BITS; i<= 1 << MAX_ADAP_D_BITS; i+=1) {
          for(j=1 << MIN_ADAP_D_BITS; j<= 1 << MAX_ADAP_D_BITS; j+=1){
            adaptiveIndexing( i, j ) ;
          }
        }
    break;

    case HV:
        
        for(i= min_size; i<= max_size; i+=1) {
          for(j=min_size; j<= max_size; j+=1){
            HV_Indexing( i, j ) ;
          }
        }

    break;
  }

 
// exit(0);
  bits_per_coordinate_w = ceil(log(image_width  / SHIFT ) / log(2.0));
  bits_per_coordinate_h = ceil(log(image_height / SHIFT ) / log(2.0));

  // for(int y=0; y < image_height / 2 - 1;y++){
  //   for(int x=0; x < image_width/2 -1; x++){
  //     printf("Contract[%d][%d] = %f\n",y,x,contract[y][x]);
  //   }
  // }

  for(k=0;k<virtual_size;k++)
  for(h=0;h<virtual_size;h++)
    qtt[k][h] = 255;

  for(k=0;k<virtual_size;k++){
    qtt[k][0] = 0;
    qtt[k][virtual_size-1] = 0;
  }

  for(k=0;k<virtual_size;k++){
    qtt[0][k] = 0;
    qtt[virtual_size-1][k] = 0;
  }

  // for(int r=0; r < virtual_size; ++r)
  //   for(int c=0; c < virtual_size; ++c)
  //     qtt[r][c] = image[r][c];

for(int r=0; r < virtual_size; ++r)
    for(int c=0; c < virtual_size; ++c)
      hv[r][c] = image[r][c];
  

  int_max_alfa = 0.5 + (MAX_ALFA )/( 8.0)* (1 << 8);
  if(int_max_alfa < 0)  int_max_alfa = 0;
  if(int_max_alfa >= (1 << 8))  int_max_alfa = (1 << 8) - 1;

  MAX_ALFA = (double) int_max_alfa / (double)(1 << 8) * ( 8.0) ;
  zeroalfa = 0;

  if ((fp = fopen(fileout, "w")) == NULL)
      fatal("\nCan't open output file");
  
  /* Header of output file */

  if(isHV){
    pack(2,(long)1,fp);
    pack(4,(long)N_BITALFA,fp);
    pack(4,(long)N_BITBETA,fp);
    // pack(4,(long)N_BITALFA1,fp);
    // pack(4,(long)N_BITALFA2,fp);
    // pack(4,(long)N_BITBETA2,fp);
  }else if(isLumInv){
    pack(2,(long)4,fp);
    pack(4,(long)N_BITALFA,fp);
    pack(4,(long)N_BITRMEAN,fp);
  }else if(isTesting){
    pack(2,(long)2,fp);
    pack(4,(long)N_BITALFA,fp);
    pack(4,(long)N_BITRMEAN,fp);
  }else if(isCovar2){
    pack(2,(long)3,fp);
    pack(4,(long)N_BITALFA,fp);
    pack(4,(long)N_BITRMEAN,fp);
  } else{
    pack(2,(long)0,fp);
    pack(4,(long)N_BITALFA,fp);
    pack(4,(long)N_BITBETA,fp);
  }
  
  

  pack(7,(long)min_size,fp);
  pack(7,(long)max_size,fp);
  pack(6,(long)SHIFT,fp);
  pack(12,(long)image_width,fp);
  pack(12,(long)image_height,fp);
  pack(8,(long) int_max_alfa,fp);
  if(isColor)
    pack(1,(long)1,fp);
  else
    pack(1,(long)0,fp);
  
  // printf("\n Image Entropy      : %f\n", entropy(image_width,image_height,0,0));
  // printf(" Image Variance     : %f\n", variance(image_width,image_height,0,0));
  // printf(" Entropy threshold  : %f\n",T_ENT);
  // printf(" Variance threshold : %f\n",T_VAR);
  // printf(" Rms threshold      : %f\n",T_RMS);
  // printf(" Color              : %s\n\n", isColor?"True":"False");

  if(isNonlinear)
    
    Nonlinear_quadtree(0,0,virtual_size,T_ENT,T_RMS,T_VAR);
  else if(isLumInv)
    LumInv_quadtree(0,0,virtual_size,T_ENT,T_RMS,T_VAR);
  else if(isTesting)
    testing_quadtree(0,0,virtual_size,T_ENT,T_RMS,T_VAR);
  else if(isCovar2)
    quadtree_2(0,0,virtual_size,T_ENT,T_RMS,T_VAR);
  else{
    if(partition_type == Quadtree)
      quadtree(0,0,virtual_size,T_ENT,T_RMS,T_VAR);
    else if(partition_type == AdaptiveQuadtree)
      traverseImage(0,0,virtual_size,virtual_size);
    else if(partition_type == HV)
      HV_traverseImage(0,0,virtual_size,virtual_size);
  }
  pack(-1,(long)0,fp);
  i = pack(-2,(long)0,fp);
  fclose(fp);
  
  t = clock() - t;
  time_taken = ((double)t)/CLOCKS_PER_SEC;
	
	
  printf("\n\n Zero_alfa_transformations   : %d\n",zero_alfa_transform);
  printf(" Number of transformations   : %d\n",transforms);
  printf(" Number of comparisons       : %ld\n",comparisons);
  printf(" Comparisons/Transformations : %f\n",(double)comparisons/transforms);
  printf(" To execute Compression tooks: %f sec  \n", time_taken);
  printf(" %d bytes written in %s\n",i,fileout);

    writeimage_pgm("hv_partition.pgm",hv,image_width,image_height); 

  if(qtree) 
    writeimage_pgm("quadtree.pgm",qtt,image_width,image_height); 


    char rows[10], cols[10];
  sprintf(rows,"%d",image_width);
  sprintf(cols,"%d",image_height);
  char *token ;
  int tok_count = 1;
  token = strtok(filein,"/");
  char *last_token;
  while(token != NULL){
   // printf("token: %s\n",token);
    last_token = token;
    token = strtok(NULL, "/");
   
  }
   // printf("token: %s\n",last_token);
  token = last_token;
  
  char image_name[strlen(token) + strlen(rows)+strlen(cols)+4];
  strcat(image_name,token);
  strcat(image_name,"(");
 
  strcat(image_name,rows);
  strcat(image_name," x ");
  strcat(image_name,cols);
  strcat(image_name,")");
  FILE *outFile = fopen("./result/as_fisher_m_hv/output.rms.csv","a+");
  fprintf(outFile,"%s\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t",image_name,SHIFT,T_RMS,time_taken,T_ENT,T_VAR,entropy(image_width,image_height,0,0),variance(image_width,image_height,0,0),COMPRESS,BPP);
  //printf("Compress: %g \t Bpp: %g\n",COMPRESS,BPP);
  // printf("%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t",image_name,SHIFT,T_RMS,time_taken,T_ENT,T_VAR,entropy(image_width,image_height,0,0),variance(image_width,image_height,0,0),COMPRESS,BPP);
  fflush(outFile);
  fclose(outFile);
  if(qtree) 
    writeimage_pgm("quadtree.pgm",qtt,image_width,image_height); 

  free(image[0]);
  free(contract[0]);
  free(flip_range[0]);
  free(range_tmp[0]);
  free(range[0]);
 
  return(0);
}



