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
  int i,k,max;
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

  Mat input_image = imread(filein, CV_LOAD_IMAGE_GRAYSCALE);
    int w = input_image.cols;
    int h = input_image.rows;
    matrix_allocate(image, w, h, PIXEL);

    for (int iii = 0; iii < h; iii++) {
        for (int jjj = 0; jjj < w; jjj++) {
            image[jjj][iii] = input_image.at<uchar>(jjj, iii);
        }
    }

  image_width = w;
    image_height = h;

  max = image_height;
  if(image_width > image_height ) 
    max = image_width; 

  virtual_size = 1 << (int) ceil(log((double) max) / log(2.0));

  matrix_allocate(contract,1+image_width / 2,1+image_height / 2,double)
  matrix_allocate(qtt,virtual_size,virtual_size,PIXEL)
  matrix_allocate(range_tmp,64,64,double)
  matrix_allocate(flip_range,64,64,double)
  matrix_allocate(range,64,64,double)

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
    case Fisher :  
      for(k=0;k<8;k++)
      for(h=0;h<3;h++)
      for(i=0;i<24;i++)
          class_fisher[k][h][i] = NULL;

      Indexing = FisherIndexing;
      Coding = FisherCoding;
      printf(" Speed-up method: Fisher\n\n");
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
    case Your_method :                               If you want to try a new 
      classification = YourMethodIndexing;           method you need just to
      coding = YourMethodCoding;                     write an indexing and a
      printf(" Speed-up method: YourMethod\n\n");    coding function and call
      break;                                         them here 
   */
  } 

  contraction(contract,image,0,0);

  for(i=(int) rint(log((double) 2) / log(2.0)); 
                     i<= (int) rint(log((double) max_size) / log(2.0)); i++) {
       Indexing((int) rint(pow(2.0,(double)i)),i);
  }

  bits_per_coordinate_w = ceil(log(image_width  / SHIFT ) / log(2.0));
  bits_per_coordinate_h = ceil(log(image_height / SHIFT ) / log(2.0));

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

  int_max_alfa = 0.5 + (MAX_ALFA )/( 8.0)* (1 << 8);
  if(int_max_alfa < 0)  int_max_alfa = 0;
  if(int_max_alfa >= (1 << 8))  int_max_alfa = (1 << 8) - 1;

  MAX_ALFA = (double) int_max_alfa / (double)(1 << 8) * ( 8.0) ;
  zeroalfa = 0;

  if ((fp = fopen(fileout, "w")) == NULL)
      fatal("\nCan't open output file");
  
  /* Header of output file */
  pack(4,(long)N_BITALFA,fp);
  pack(4,(long)N_BITBETA,fp);
  pack(7,(long)min_size,fp);
  pack(7,(long)max_size,fp);
  pack(6,(long)SHIFT,fp);
  pack(12,(long)image_width,fp);
  pack(12,(long)image_height,fp);
  pack(8,(long) int_max_alfa,fp);

  printf("\n Image Entropy      : %f\n", entropy(image_width,image_height,0,0));
  printf(" Image Variance     : %f\n", variance(image_width,image_height,0,0));
  printf(" Entropy threshold  : %f\n",T_ENT);
  printf(" Variance threshold : %f\n",T_VAR);
  printf(" Rms threshold      : %f\n\n",T_RMS);

  quadtree(0,0,virtual_size,T_ENT,T_RMS,T_VAR);

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

  if(qtree) 
    writeimage_pgm("quadtree.pgm",qtt,image_width,image_height); 

  free(image[0]);
  free(contract[0]);
  free(flip_range[0]);
  free(range_tmp[0]);
  free(range[0]);
 
  return(0);
}



