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
#include <unistd.h>
#include <time.h>
#include "def.h"
#define EXTERN
#include "globals.h"
#include "triangle.h"
#include "prot.h"
//#include <vector>


#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include "GMRsaliency.h"



int main(int argc, char **argv)
{
  int int_max_alfa,step;
  int pipe_disp[2], pid;
  FILE *output;

  printf("%s",triangle_4);

  getopt_dec(argc, argv);

  if ((input = fopen(filein, "r")) == NULL)
      fatal("\n Can't open input file");

  unpack(-2,input);    /*Initialize unpack */
  int type    = (int)unpack(2,input);
  
  switch(type){
    case 0:

      N_BITALFA = (int)unpack(4,input);
      N_BITBETA = (int)unpack(4,input);
      break;
    case 1:
      printf("Fisher Coding with HV partition");
      isHV = 1;
      N_BITALFA = (int)unpack(4,input);
      N_BITBETA = (int)unpack(4,input);
      break;
      // isNonlinear = 1;
      // printf("Nonlinear Fisher Coding\n");
      // N_BITALFA1 =  (int)unpack(4,input);
      // N_BITALFA2 =  (int)unpack(4,input);
      // N_BITBETA2 =  (int)unpack(4,input);
      // break;
    case 4:
      isLumInv = 1;
      printf("Luminance Invariant Fisher coding\n");
      N_BITALFA = (int)unpack(4,input);
      N_BITRMEAN = (int)unpack(4,input);
      break;
    case 2:
      isTesting = 1;
      printf("testing new method with Fisher coding\n");
      N_BITALFA = (int)unpack(4,input);
      N_BITRMEAN = (int)unpack(4,input);
      break;
    case 3:
      isCovar2 = 1;
      printf("Cov class, adaptive search and adaptive range quadtree\n");
      N_BITALFA = (int)unpack(4,input);
      N_BITRMEAN = (int)unpack(4,input);
      break;
  }
  
  min_size       = (int)unpack(7,input);
  max_size       = (int)unpack(7,input);
  SHIFT          = (int)unpack(6,input);
  image_width   = (int)unpack(12,input);
  image_height   = (int)unpack(12,input);
  int_max_alfa   = (int)unpack(8,input);
  isColor        = (int)unpack(1,input);

  printf("Color: %s\n", isColor?"True":"False");

  bits_per_coordinate_w = ceil(log(image_width  / SHIFT ) / log(2.0));
  bits_per_coordinate_h = ceil(log(image_height / SHIFT ) / log(2.0));

  zeroalfa = 0;
  MAX_ALFA = (double) int_max_alfa / (double)(1 << 8) * ( 8.0) ;

  max = image_height;
  min = image_width;
  if(image_width > image_height ) {
    min = image_height;
    max = image_width;
  }
  
  virtual_size = 1 << (int) ceil(log((double) max) / log(2.0));

  trans = &fractal_code; 
  printf("\n Reading %s ... ",filein);
  fflush(stdout);

  printf(" Original image size: %dx%d\n",image_width,image_height);

  image_width  = (int) rint((zoom * image_width));
  image_height = (int) rint((zoom * image_height));

  if(zoom != 1.0) {
    printf(" Zooming image to   : %dx%d\n",image_width,image_height);
    fflush(stdout);
  } 

  matrix_allocate(image,2+image_width,2+image_height,PIXEL)
  matrix_allocate(image1,2+image_width,2+image_height,PIXEL)
  if(isColor){
    matrix_allocate(image_uch,2+image_width,2+image_height,PIXEL)
    matrix_allocate(image_vch,2+image_width,2+image_height,PIXEL)
  }

  Mat img = imread("lenna.pgm", CV_LOAD_IMAGE_GRAYSCALE);
  GMRsaliency GMRsal;
  Mat sal=GMRsal.GetSal(img);

  imwrite("sal_image.png",sal*255);

  if(isNonlinear){
    read_transformations_nonlinear(0,0,virtual_size);
  }else if(isLumInv){ 
    read_transformations_LumInv(0,0,virtual_size);
  }else if(isTesting){
    read_transformations_testing(0,0,virtual_size);
  }else if(isCovar2){
    read_transformations_2(0,0,virtual_size);
  }else if(iterDec2){
    printf(" Decoding Using proposed initial image\n");
    read_transformations_new_init_image(0,0,virtual_size);
  }else if(isHV){
    printf(" Reading HV transforms...\n");
    HV_traverseImage_2(0,0,virtual_size,virtual_size);    
  }else{
    // traverseImage_2(0,0,virtual_size,virtual_size);
    read_transformations(0,0,virtual_size);

  }
  
  printf(" Tranforms = %d\n",transforms);
  printf(" done\n");
  fflush(stdout);

  // printf(" Original image size: %dx%d\n",image_width,image_height);

  // image_width  = (int) rint((zoom * image_width));
  // image_height = (int) rint((zoom * image_height));

  // if(zoom != 1.0) {
  //   printf(" Zooming image to   : %dx%d\n",image_width,image_height);
  //   fflush(stdout);
  // } 

  // matrix_allocate(image,2+image_width,2+image_height,PIXEL)
  // matrix_allocate(image1,2+image_width,2+image_height,PIXEL)
  // if(isColor){
  //   matrix_allocate(image_uch,2+image_width,2+image_height,PIXEL)
  //   matrix_allocate(image_vch,2+image_width,2+image_height,PIXEL)
  // }
  if(piramidal) {
      min *= zoom;
      step = SHIFT * floor(zoom);
      if(step == 0) step = 1;
      lev = 0;
      while(1){
        if(min < 200 || (step & 1))
            break;
        min  >>= 1;
        step >>= 1;
        lev++;
      }
     printf("\n %d level piramid\n",lev);
     
     if(isNonlinear){
      iterative_decoding_nonlinear(lev,iterations,zoom);    /* Decode at low resolution */ 
     }else if(isLumInv){
      iterative_decoding_LumInv(lev,iterations,zoom);
     }else if(isTesting){
      iterative_decoding_testing(lev,iterations,zoom);
     }else if(iterDec2){
      iterative_decoding_new_init_image(lev,iterations,zoom);
     }else{
      iterative_decoding(lev,iterations,zoom);
     }


     if(isNonlinear){ 
       piramidal_decoding_nonlinear(lev); /* Increase resolution      */
     }else if(isLumInv){
       piramidal_decoding_LumInv(lev);                     
     }else if(isTesting){
       piramidal_decoding_testing(lev);
     } else{
       piramidal_decoding(lev);                     
     }
       
     
     if(quality){
        if(isNonlinear)
          iterative_decoding_nonlinear(0,2,1.0); 
        else if (isLumInv)
          iterative_decoding_LumInv(0,2,1.0);
        else if (isTesting)
          iterative_decoding_testing(0,2,1.0);
        else if (iterDec2)
          iterative_decoding_new_init_image(0,2,1.0);
         else 
          iterative_decoding(0,2,1.0);     

     }    
  
  } else{
      if(isNonlinear)
        iterative_decoding_nonlinear(0,iterations,zoom);
      else if(isLumInv)
        iterative_decoding_LumInv(0,iterations,zoom); 
      else if(isTesting)
        iterative_decoding_testing(0,iterations,zoom);
      else if(iterDec2)
        iterative_decoding_new_init_image(0,iterations,zoom);
      else
        iterative_decoding(0,iterations,zoom); 

    }
  

  if(postproc) 
     smooth_image(); 

  if( isColor ){ // Conver to RGB
    std::vector<Mat> yuvChannels;
    Mat ych(image_height,image_width, CV_8U);
    Mat uch(image_height,image_width, CV_8U);
    Mat vch(image_height,image_width, CV_8U);
    yuvChannels.push_back(ych);
    yuvChannels.push_back(uch);
    yuvChannels.push_back(vch);

    for (int iii = 0; iii < image_width; iii++) {
      for (int jjj = 0; jjj < image_height; jjj++) {
        yuvChannels[0].at<uchar>(jjj,iii) = image[jjj][iii];
        yuvChannels[1].at<uchar>(jjj,iii) = image_uch[jjj][iii];
        yuvChannels[2].at<uchar>(jjj,iii) = image_vch[jjj][iii];
      }
    }
    merge(yuvChannels, ych);
    cvtColor(ych, ych, CV_YUV2BGR);
    printf("Writing... output.dec.ppm");
    imwrite(fileout, ych);
    
  }

  if(display)  {
     if ( pipe(pipe_disp)) 
        fatal("\n Pipe creation error");

     if ((pid=fork()) < 0)
        fatal("\n Fork error");

     if (pid==0) {     /* Child */
         dup2(pipe_disp[0],fileno(stdin));
         close(pipe_disp[0]); 
         close(pipe_disp[1]); 
	 execlp("xv","xv","-",(char *) 0);
         fatal("\n Exec error xv not started");
     }

     printf("\n");
     close(pipe_disp[0]); 
     output = fdopen(pipe_disp[1], "w");
     writeimage_pipe(output, image,image_width,image_height);
  } 
  else 
    if(raw_format)  {}
       // writeimage_raw(fileout, image,image_width,image_height);
    else{}
      // writeimage_pgm(fileout, image,image_width,image_height);
    
    
     
  free(image[0]);
  free(image1[0]);

  return 0;
}




void zooming(double scalefactor)
{
  trans = &fractal_code;
  while (trans->next != NULL) {
     trans = trans->next;

     trans->rrx  *= scalefactor;
     trans->rry  *= scalefactor;
     trans->rx   *= scalefactor;
     trans->ry   *= scalefactor;
     trans->dx   *= scalefactor;
     trans->dy   *= scalefactor;
  }
}

void HV_traverseImage_2(int atx, int aty, int x_size, int y_size)
{
  
  int s_size;
  int s_log;

  int max_bits = (int) log2(max_size);
  s_log = (int) log2(min_2(x_size,y_size));
  s_size = 1 << s_log;
  // printf("s_size = %d\n",s_size);

  // printf("s_log = %d\n",s_log);
  if(s_log > max_bits){
    HV_traverseImage_2(atx,          aty,          s_size/2, s_size/2);
    HV_traverseImage_2(atx+s_size/2, aty,          s_size/2, s_size/2);
    HV_traverseImage_2(atx,          aty+s_size/2, s_size/2, s_size/2);
    HV_traverseImage_2(atx+s_size/2, aty+s_size/2, s_size/2, s_size/2);
  }
  else{
    HV_decompressRange(atx,aty, x_size, y_size);
  }

  if(x_size > s_size)
    HV_traverseImage_2(atx+s_size,aty,x_size - s_size, y_size);
  
  if(y_size > s_size)
    HV_traverseImage_2(atx, aty + s_size,s_size, y_size - s_size);
}

void HV_decompressRange(int atx,int aty,int x_size,int y_size)
{ 
  int qalfa,qbeta;
  double alfa,beta, um,vm;
  int ddx, ddy;


  // if(atx >= image_height  || aty >= image_width )
  //     return;
  
  if(x_size==0||y_size==0)
    return;

  if(min_2(x_size,y_size) == min_size || (int)unpack(1,input) == 1){
    transforms++;
    trans->next = (struct t_node *) malloc(sizeof(struct t_node ));
    trans       = trans->next; 
    trans->next = NULL;
    qalfa       = (int)unpack(N_BITALFA,  input);
    qbeta       = (int)unpack(N_BITBETA,  input);
    // if(isColor){
    //   um = (int)unpack(8,input);
    //   vm = (int)unpack(8,input);
    // }

    /* Compute alfa from the quantized value */
    alfa = (double) qalfa / (double)(1 << N_BITALFA) * (MAX_ALFA) ;
    
    /* Compute beta from the quantized value */
    beta = (double) qbeta/(double)((1 << N_BITBETA)-1)* ((1.0+fabs(alfa)) * 255);
    if (alfa > 0.0) beta  -= alfa * 255;
       
    trans->alfa = alfa;
    trans->beta = beta;
    // if(isColor){
    //   trans->um = um;
    //   trans->vm = vm;
    // }
    if(qalfa != zeroalfa) {      
        trans-> sym_op = (int)unpack(3, input);
        ddx = (int)unpack(bits_per_coordinate_h,input);
        ddy = (int)unpack(bits_per_coordinate_w,input);
        trans->dx = SHIFT * ddx;
        trans->dy = SHIFT * ddy;
       // printf("%d %d\n", ddx, ddy);
    } else {
        trans-> sym_op = 0;
        trans-> dx  = 0;
        trans-> dy = 0;
    }
    trans->rx = atx;
    trans->ry = aty;
    trans->x_size = x_size;
    trans->y_size = y_size;
    trans->rrx = atx + x_size;
    trans->rry = aty + y_size;

  } else {
    int l1,l2,l3;

    l1 = (int)unpack(1,input);
    l2 = (int)unpack(3,input);

    if(l1){
      HV_decompressRange(atx, aty,      x_size, l2);
      HV_decompressRange(atx, aty + l2, x_size, y_size-l2);
    }else{
      HV_decompressRange(atx,      aty, l2,        y_size);
      HV_decompressRange(atx + l2, aty, x_size-l2, y_size);
    }
    
  } 
}

void read_transformations_testing(int atx,int aty,int size)
{ 
  int qalfa,qrmean;
  double alfa,rmean,beta, um,vm;
  int ddx, ddy;


  if(atx >= image_height  || aty >= image_width )
      return;

  if (size > max_size || atx+size > image_height || aty+size > image_width){
      read_transformations_testing(atx,aty,size/2);
      read_transformations_testing(atx+size/2,aty,size/2);
      read_transformations_testing(atx,aty+size/2,size/2);
      read_transformations_testing(atx+size/2,aty+size/2,size/2);
      return;
  }

  if (size > min_size && unpack(1,input)) {
      /* A 1 means we subdivided.. so quadtree */
      read_transformations_testing(atx,aty,size/2);
      read_transformations_testing(atx+size/2,aty,size/2);
      read_transformations_testing(atx,aty+size/2,size/2);
      read_transformations_testing(atx+size/2,aty+size/2,size/2);
  } else {
      /* Read the trasformation */  
      transforms++;
      trans->next = (struct t_node *) malloc(sizeof(struct t_node ));
      trans       = trans->next; 
      trans->next = NULL;
      qalfa       = (int)unpack(N_BITALFA,  input);
      qrmean       = (int)unpack(N_BITRMEAN,  input);
      if(isColor){
        um = (int)unpack(8,input);
        vm = (int)unpack(8,input);
      }

      /* Compute alfa from the quantized value */
      alfa = (double) qalfa / (double)(1 << N_BITALFA) * (MAX_ALFA) ;
      
      /* Compute rmean from the quantized value */
     rmean = (double)qrmean/(double)((1 << N_BITRMEAN)-1)*((1.0+fabs(alfa))*255);
        if (alfa > 0.0) rmean  -= alfa * 255;
         
      trans->alfa = alfa;
      trans->rmean = rmean;
      if(isColor){
        trans->um = um;
        trans->vm = vm;
      }
      if(qalfa != zeroalfa) {      
          trans-> sym_op = (int)unpack(3, input);
          ddx = (int)unpack(bits_per_coordinate_h,input);
          ddy = (int)unpack(bits_per_coordinate_w,input);
          trans->dx = SHIFT * ddx;
          trans->dy = SHIFT * ddy;
         // printf("%d %d\n", ddx, ddy);
      } else {
          trans-> sym_op = 0;
          trans-> dx  = 0;
          trans-> dy = 0;
      }
      trans->rx = atx;
      trans->ry = aty;
      trans->size = size;
      trans->rrx = atx + size;
      trans->rry = aty + size;
 
  }
}
void read_transformations_LumInv(int atx,int aty,int size)
{ 
  int qalfa,qrmean;
  double alfa,rmean,beta, um,vm;
  int ddx, ddy;


  if(atx >= image_height  || aty >= image_width )
      return;

  if (size > max_size || atx+size > image_height || aty+size > image_width){
      read_transformations_LumInv(atx,aty,size/2);
      read_transformations_LumInv(atx+size/2,aty,size/2);
      read_transformations_LumInv(atx,aty+size/2,size/2);
      read_transformations_LumInv(atx+size/2,aty+size/2,size/2);
      return;
  }

  if (size > min_size && unpack(1,input)) {
      /* A 1 means we subdivided.. so quadtree */
      read_transformations_LumInv(atx,aty,size/2);
      read_transformations_LumInv(atx+size/2,aty,size/2);
      read_transformations_LumInv(atx,aty+size/2,size/2);
      read_transformations_LumInv(atx+size/2,aty+size/2,size/2);
  } else {
      /* Read the trasformation */  
      transforms++;
      trans->next = (struct t_node *) malloc(sizeof(struct t_node ));
      trans       = trans->next; 
      trans->next = NULL;
      qalfa       = (int)unpack(N_BITALFA,  input);
      qrmean       = (int)unpack(N_BITRMEAN,  input);
      if(isColor){
        um = (int)unpack(8,input);
        vm = (int)unpack(8,input);
      }

      /* Compute alfa from the quantized value */
      alfa = (double) qalfa / (double)(1 << N_BITALFA) * (MAX_ALFA) ;
      
      /* Compute rmean from the quantized value */
     rmean = (double)qrmean/(double)((1 << N_BITRMEAN)-1)*((1.0+fabs(alfa))*255);
        if (alfa > 0.0) beta  -= alfa * 255;
         
      trans->alfa = alfa;
      trans->rmean = rmean;
      if(isColor){
        trans->um = um;
        trans->vm = vm;
      }
      if(qalfa != zeroalfa) {      
          trans-> sym_op = (int)unpack(3, input);
          ddx = (int)unpack(bits_per_coordinate_h,input);
          ddy = (int)unpack(bits_per_coordinate_w,input);
          trans->dx = SHIFT * ddx;
          trans->dy = SHIFT * ddy;
         // printf("%d %d\n", ddx, ddy);
      } else {
          trans-> sym_op = 0;
          trans-> dx  = 0;
          trans-> dy = 0;
      }
      trans->rx = atx;
      trans->ry = aty;
      trans->size = size;
      trans->rrx = atx + size;
      trans->rry = aty + size;
 
  }
}

void read_transformations_nonlinear(int atx,int aty,int size)
{ 
  int qalfa,qalfa1,qalfa2,qbeta;
  double alfa,alfa1,alfa2,beta, um,vm;
  int ddx, ddy;

  if(atx >= image_height  || aty >= image_width )
      return;

  if (size > max_size || atx+size > image_height || aty+size > image_width){
      read_transformations_nonlinear(atx,aty,size/2);
      read_transformations_nonlinear(atx+size/2,aty,size/2);
      read_transformations_nonlinear(atx,aty+size/2,size/2);
      read_transformations_nonlinear(atx+size/2,aty+size/2,size/2);
      return;
  }

  if (size > min_size && unpack(1,input)) {
      /* A 1 means we subdivided.. so quadtree */
      read_transformations_nonlinear(atx,aty,size/2);
      read_transformations_nonlinear(atx+size/2,aty,size/2);
      read_transformations_nonlinear(atx,aty+size/2,size/2);
      read_transformations_nonlinear(atx+size/2,aty+size/2,size/2);
  } else {
      /* Read the trasformation */  
      transforms++;
      trans->next = (struct t_node *) malloc(sizeof(struct t_node ));
      trans       = trans->next; 
      trans->next = NULL;
      qalfa1       = (int)unpack(N_BITALFA1,  input);
      qalfa2       = (int)unpack(N_BITALFA2,  input);
      qbeta       = (int)unpack(N_BITBETA2,  input);
      if(isColor){
        um = (int)unpack(8,input);
        vm = (int)unpack(8,input);
      }

      /* Compute alfa from the quantized value */
<<<<<<< HEAD
       // alfa1 = (qalfa1 - 15) * 0.1;
       // alfa2 = (qalfa2 - 15) * 0.001;
      alfa1 = (double) qalfa1 / (double)(1 << N_BITALFA1) * ( MAX_ALFA1) ;
      alfa2 = (double) qalfa2 / (double)(1 << N_BITALFA2) * ( MAX_ALFA2) ;
=======
       alfa1 = (qalfa1 - 15) * 0.1;
       alfa2 = (qalfa2 - 15) * 0.001;
      
>>>>>>> 4d466d691a1cbf5bc45d429a910f079a7d690d53

      
      /* Compute beta from the quantized value */
      // beta = (qbeta - 128);
       beta = (double)qbeta/(double)((1 << N_BITBETA2)-1)*((1.0+fabs(alfa1))*255);
      if (alfa1 > 0.0) beta  -= alfa1 * 255;
      // alfa2 = (double) qalfa2 / (double)(1 << N_BITALFA2) * ( MAX_ALFA) ;
      // beta = (double)qbeta/(double)((1 << N_BITBETA2)-1)*((1.0+fabs(alfa1))*255);
      // if (alfa1 > 0.0) beta  -= alfa1 * 255;
         // printf("Alfa1 = %f \t Alfa2 = %f \n",alfa1,alfa2);
      trans->alfa1 = alfa1;
      trans->alfa2 = alfa2;
      trans->beta = beta;
      if(isColor){
        trans->um = um;
        trans->vm = vm;
      }
     if(qalfa1 != zeroalfa || qalfa2 != zeroalfa) {      
          trans-> sym_op = (int)unpack(3, input);
          ddx = (int)unpack(bits_per_coordinate_h,input);
          ddy = (int)unpack(bits_per_coordinate_w,input);
          trans->dx = SHIFT * ddx;
          trans->dy = SHIFT * ddy;
         // printf("%d %d\n", ddx, ddy);
     } else {
          trans-> sym_op = 0;
          trans-> dx  = 0;
          trans-> dy = 0;
      }
      trans->rx = atx;
      trans->ry = aty;
      trans->size = size;
      trans->rrx = atx + size;
      trans->rry = aty + size;
 
  }
}

void read_transformations_2(int atx,int aty,int size)
{ 
  int qalfa,qbeta;
  double alfa,beta, um,vm;
  int ddx, ddy;

  int s_log, s_size;
  s_log = bitlength(size) - 1;
  s_size = 1 << s_log;


  if(atx >= image_height  || aty >= image_width )
      return;

  if (size > max_size || atx+size > image_height || aty+size > image_width){
      read_transformations_2(atx,aty,size/2);
      read_transformations_2(atx+size/2,aty,size/2);
      read_transformations_2(atx,aty+size/2,size/2);
      read_transformations_2(atx+size/2,aty+size/2,size/2);
      return;
  }

  if (size > min_size && unpack(1,input)) {
      /* A 1 means we subdivided.. so quadtree */
      read_transformations_2(atx,aty,size/2);
      read_transformations_2(atx+size/2,aty,size/2);
      read_transformations_2(atx,aty+size/2,size/2);
      read_transformations_2(atx+size/2,aty+size/2,size/2);
  }else if(s_log > MAX_BITS && unpack(1,input)){
      read_transformations_2(atx,aty,size/2);
      read_transformations_2(atx+size/2,aty,size/2);
      read_transformations_2(atx,aty+size/2,size/2);
      read_transformations_2(atx+size/2,aty+size/2,size/2);
  } else {
      /* Read the trasformation */  
      transforms++;
      trans->next = (struct t_node *) malloc(sizeof(struct t_node ));
      trans       = trans->next; 
      trans->next = NULL;
      qalfa       = (int)unpack(N_BITALFA,  input);
      qbeta       = (int)unpack(N_BITBETA,  input);
      if(isColor){
        um = (int)unpack(8,input);
        vm = (int)unpack(8,input);
      }

      /* Compute alfa from the quantized value */
      alfa = (double) qalfa / (double)(1 << N_BITALFA) * (MAX_ALFA) ;
      
      /* Compute beta from the quantized value */
      beta = (double) qbeta/(double)((1 << N_BITBETA)-1)* ((1.0+fabs(alfa)) * 255);
      if (alfa > 0.0) beta  -= alfa * 255;
         
      trans->alfa = alfa;
      trans->beta = beta;
      if(isColor){
        trans->um = um;
        trans->vm = vm;
      }
      if(qalfa != zeroalfa) {      
          trans-> sym_op = (int)unpack(3, input);
          ddx = (int)unpack(bits_per_coordinate_h,input);
          ddy = (int)unpack(bits_per_coordinate_w,input);
          trans->dx = SHIFT * ddx;
          trans->dy = SHIFT * ddy;
         // printf("%d %d\n", ddx, ddy);
      } else {
          trans-> sym_op = 0;
          trans-> dx  = 0;
          trans-> dy = 0;
      }
      trans->rx = atx;
      trans->ry = aty;
      trans->size = size;
      trans->rrx = atx + size;
      trans->rry = aty + size;
 
  }
}

void traverseImage_2(int atx, int aty, int x_size, int y_size)
{
  
  int s_size;
  int s_log;

  s_log = (int) log2(min_2(x_size,y_size));
  s_size = 1 << s_log;
  // printf("s_size = %d\n",s_size);

  // printf("s_log = %d\n",s_log);
  if(s_log > MAX_ADAP_R_BITS){
    traverseImage_2(atx,          aty,          s_size/2, s_size/2);
    traverseImage_2(atx+s_size/2, aty,          s_size/2, s_size/2);
    traverseImage_2(atx,          aty+s_size/2, s_size/2, s_size/2);
    traverseImage_2(atx+s_size/2, aty+s_size/2, s_size/2, s_size/2);
  }
  else{
    decompressRange(atx,aty, x_size, y_size);
  }

  if(x_size > s_size)
    traverseImage_2(atx+s_size,aty,x_size - s_size, y_size);
  
  if(y_size > s_size)
    traverseImage_2(atx, aty + s_size,s_size, y_size - s_size);
}

void decompressRange(int atx,int aty,int x_size,int y_size)
{ 
  int qalfa,qbeta;
  double alfa,beta, um,vm;
  int ddx, ddy;


  // if(atx >= image_height  || aty >= image_width )
  //     return;
  
  if(x_size==0||y_size==0)
    return;

  if(max_2(x_size,y_size) == 1 << MIN_ADAP_R_BITS || (int)unpack(1,input) == 1){
    transforms++;
    trans->next = (struct t_node *) malloc(sizeof(struct t_node ));
    trans       = trans->next; 
    trans->next = NULL;
    qalfa       = (int)unpack(N_BITALFA,  input);
    qbeta       = (int)unpack(N_BITBETA,  input);
    // if(isColor){
    //   um = (int)unpack(8,input);
    //   vm = (int)unpack(8,input);
    // }

    /* Compute alfa from the quantized value */
    alfa = (double) qalfa / (double)(1 << N_BITALFA) * (MAX_ALFA) ;
    
    /* Compute beta from the quantized value */
    beta = (double) qbeta/(double)((1 << N_BITBETA)-1)* ((1.0+fabs(alfa)) * 255);
    if (alfa > 0.0) beta  -= alfa * 255;
       
    trans->alfa = alfa;
    trans->beta = beta;
    // if(isColor){
    //   trans->um = um;
    //   trans->vm = vm;
    // }
    // if(qalfa != zeroalfa) {      
        trans-> sym_op = (int)unpack(3, input);
        ddx = (int)unpack(bits_per_coordinate_h,input);
        ddy = (int)unpack(bits_per_coordinate_w,input);
        trans->dx = SHIFT * ddx;
        trans->dy = SHIFT * ddy;
    //    // printf("%d %d\n", ddx, ddy);
    // } else {
    //     trans-> sym_op = 0;
    //     trans-> dx  = 0;
    //     trans-> dy = 0;
    // }
    trans->rx = atx;
    trans->ry = aty;
    trans->x_size = x_size;
    trans->y_size = y_size;
    trans->rrx = atx + x_size;
    trans->rry = aty + y_size;

  } else {
    int l1,l2;
    //  printf("%s\n","in adaptive part");
    l1 = (int)unpack(3,input);
    l2 = (int)unpack(3,input);
     
     // l1 = 1 << (int)log2(l1);   
     // l2 = 1 << (int)log2(l2);   
     // l3= 1 << (int)log2(l3);   
    
    decompressRange(atx,      aty,      l2,      l1);
    decompressRange(atx+l2,  aty,      x_size-l2, l1);
    decompressRange(atx,      aty+l1,  l2,      y_size-l1);
    decompressRange(atx+l2,  aty+l1,   x_size-l2, y_size-l1);
  } 
}

void read_transformations_new_init_image(int atx,int aty,int size)
{ 
  int qalfa,qbeta;
  double alfa,beta, um,vm;
  int ddx, ddy;


  if(atx >= image_height  || aty >= image_width )
      return;

  if (size > max_size || atx+size > image_height || aty+size > image_width){
      read_transformations_new_init_image(atx,aty,size/2);
      read_transformations_new_init_image(atx+size/2,aty,size/2);
      read_transformations_new_init_image(atx,aty+size/2,size/2);
      read_transformations_new_init_image(atx+size/2,aty+size/2,size/2);
      return;
  }

  if (size > min_size && unpack(1,input)) {
      /* A 1 means we subdivided.. so quadtree */
      read_transformations_new_init_image(atx,aty,size/2);
      read_transformations_new_init_image(atx+size/2,aty,size/2);
      read_transformations_new_init_image(atx,aty+size/2,size/2);
      read_transformations_new_init_image(atx+size/2,aty+size/2,size/2);
  } else {
      /* Read the trasformation */  
      transforms++;
      trans->next = (struct t_node *) malloc(sizeof(struct t_node ));
      trans       = trans->next; 
      trans->next = NULL;
      qalfa       = (int)unpack(N_BITALFA,  input);
      qbeta       = (int)unpack(N_BITBETA,  input);
      if(isColor){
        um = (int)unpack(8,input);
        vm = (int)unpack(8,input);
      }

      /* Compute alfa from the quantized value */
      alfa = (double) qalfa / (double)(1 << N_BITALFA) * (MAX_ALFA) ;
      
      /* Compute beta from the quantized value */
      beta = (double) qbeta/(double)((1 << N_BITBETA)-1)* ((1.0+fabs(alfa)) * 255);
      if (alfa > 0.0) beta  -= alfa * 255;
         
      trans->alfa = alfa;
      trans->beta = beta;
      if(isColor){
        trans->um = um;
        trans->vm = vm;
      }
      if(qalfa != zeroalfa) {      
          trans-> sym_op = (int)unpack(3, input);
          ddx = (int)unpack(bits_per_coordinate_h,input);
          ddy = (int)unpack(bits_per_coordinate_w,input);
          trans->dx = SHIFT * ddx;
          trans->dy = SHIFT * ddy;
         // printf("%d %d\n", ddx, ddy);
      } else {
          trans-> sym_op = 0;
          trans-> dx  = 0;
          trans-> dy = 0;
      }
      trans->rx = atx;
      trans->ry = aty;
      trans->size = size;
      trans->rrx = atx + size;
      trans->rry = aty + size;
      int rx, ry,rrx,rry;
      rx = (int)rint(trans->rx);
      ry = (int)rint(trans->ry);
      rrx = (int)rint(trans->rrx);
      rry= (int)rint(trans->rry);
      int val = rint(beta/(1-alfa));
      // printf("image = %f\n",beta/(1-alfa));
      for(int i=rx; i <= rrx; i++){
        for(int j=ry; j <= rry; j++){
          if(val < 0) image[i][j] = 0;
          if(val > 255) image[i][j] = 255;
          else image[i][j] = val;
 
        }
      }
      // image[dx][dy] = (unsigned char)rint(abs(beta/(1-alfa)));
      // printf("image = %u\n",(unsigned int)image[dx][dy]);
      
  }
}


void read_transformations(int atx,int aty,int size)
{ 
  int qalfa,qbeta;
  double alfa,beta, um,vm;
  int ddx, ddy;


  if(atx >= image_height  || aty >= image_width )
      return;

  if (size > max_size || atx+size > image_height || aty+size > image_width){
      read_transformations(atx,aty,size/2);
      read_transformations(atx+size/2,aty,size/2);
      read_transformations(atx,aty+size/2,size/2);
      read_transformations(atx+size/2,aty+size/2,size/2);
      return;
  }

  if (size > min_size && unpack(1,input)) {
      /* A 1 means we subdivided.. so quadtree */
      read_transformations(atx,aty,size/2);
      read_transformations(atx+size/2,aty,size/2);
      read_transformations(atx,aty+size/2,size/2);
      read_transformations(atx+size/2,aty+size/2,size/2);
  } else {
      /* Read the trasformation */  
      transforms++;
      trans->next = (struct t_node *) malloc(sizeof(struct t_node ));
      trans       = trans->next; 
      trans->next = NULL;
      qalfa       = (int)unpack(N_BITALFA,  input);
      qbeta       = (int)unpack(N_BITBETA,  input);
      if(isColor){
        um = (int)unpack(8,input);
        vm = (int)unpack(8,input);
      }

      /* for progressive decoding */
      int pqalfa = 0,pqbeta = 0;
      for(int b=0; b < N_BITS; ++b){
        if(N_BITS < N_BITALFA)
          pqalfa |= qalfa & (1 << (N_BITALFA-b));
        else
          pqalfa = qalfa;

        if(N_BITS < N_BITBETA)
          pqbeta |= qbeta & (1 << (N_BITBETA-b));
        else
          pqbeta = qbeta;
      }

      /* Compute alfa from the quantized value */
      alfa = (double) pqalfa / (double)(1 << N_BITALFA) * (MAX_ALFA) ;
      
      /* Compute beta from the quantized value */
      beta = (double) pqbeta/(double)((1 << N_BITBETA)-1)* ((1.0+fabs(alfa)) * 255);
      if (alfa > 0.0) beta  -= alfa * 255;
         
      trans->alfa = alfa;
      trans->beta = beta;
      if(isColor){
        trans->um = um;
        trans->vm = vm;
      }
      if(qalfa != zeroalfa) {      
          trans-> sym_op = (int)unpack(3, input);
          ddx = (int)unpack(bits_per_coordinate_h,input);
          ddy = (int)unpack(bits_per_coordinate_w,input);
          trans->dx = SHIFT * ddx;
          trans->dy = SHIFT * ddy;
         // printf("%d %d\n", ddx, ddy);
      } else {
          trans-> sym_op = 0;
          trans-> dx  = 0;
          trans-> dy = 0;
      }
      trans->rx = atx;
      trans->ry = aty;
      trans->size = size;
      trans->rrx = atx + size;
      trans->rry = aty + size;
      int dx, dy;
      dx = (int)rint(trans->dx);
      dy = (int)rint(trans->dy);
      
  }
}

void iterative_decoding_testing(int level,int n_iter,double zoo)
{
  int rx,ry,rrx,rry,dx,dy;
  register int i,j,ii,jj,s;
  register PIXEL **imag,**imag1;
  double pixel;
  double z_factor;
  int width,height;
  static int first_time = 0;


  printf("\n");
  z_factor = zoo / (double) (1 << level);
  zooming(z_factor);
  width  = (int) rint(image_width  * z_factor / zoo); 
  height = (int) rint(image_height * z_factor / zoo); 

  if(first_time++ == 0)
     for(i=0;i< height;i++)
     for(j=0;j< width ;j++) 
        image[i][j] = 128;

  for(s=0; s < n_iter ; s++) {
     imag = image;
     imag1 = image1;

     if(level > 0)
         printf(" Decoding at low resolution (%dx%d) %d\r",width,height,s);
     else
         printf(" Iterative decoding... %d\r",s);

     fflush(stdout);
     trans = &fractal_code;
     while(trans->next != NULL)  {
        trans = trans->next;
        rx=(int)rint(trans->rx);
        ry=(int)rint(trans->ry);
        dx=(int)rint(trans->dx);
        dy=(int)rint(trans->dy);
        rrx=(int)rint(trans->rrx);
        rry=(int)rint(trans->rry);

        switch(trans->sym_op) {     
         case IDENTITY   : 
            for(i=rx,ii=dx;i< rrx;i++,ii+=2)
            for(j=ry,jj=dy;j< rry;j++,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               trans->beta = trans->rmean - trans->alfa * pixel; 
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
               // printf("%d\n",imag1[i][j]);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case R_ROTATE90 :  
            for(j=rry-1,ii=dx;j>= ry;j--,ii+=2)
            for(i=rx,jj=dy;i< rrx;i++,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               trans->beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
               // printf("%d\n",imag1[i][j]);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }

            }
            break;
         case L_ROTATE90 :  
      for(j=ry,ii=dx;j< rry;j++,ii+=2) 
            for(i=rrx-1,jj=dy;i>= rx;i--,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               trans->beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
               // printf("%d\n",imag1[i][j]);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }

            }
            break;
         case ROTATE180 :  
      for(i=rrx-1,ii=dx;i>= rx;i--,ii+=2) 
            for(j=rry-1,jj=dy;j>= ry;j--,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               trans->beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
               // printf("%d\n",imag1[i][j]);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case R_VERTICAL :  
      for(i=rx,ii=dx;i< rrx;i++,ii+=2) 
            for(j=rry-1,jj=dy;j>= ry;j--,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               trans->beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
               // printf("%d\n",imag1[i][j]);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case R_HORIZONTAL: 
      for(i=rrx-1,ii=dx;i>= rx;i--,ii+=2)
            for(j=ry,jj=dy;j< rry;j++,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               trans->beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
               // printf("%d\n",imag1[i][j]);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case F_DIAGONAL :  
      for(j=ry,ii=dx;j< rry;j++,ii+=2)
            for(i=rx,jj=dy;i< rrx;i++,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               trans->beta = trans->rmean - trans->alfa * pixel;
               // printf("%f\n",trans->beta);
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case S_DIAGONAL:   
      for(j=rry-1,ii=dx;j>= ry;j--,ii+=2) 
      for(i=rrx-1,jj=dy;i>= rx;i--,jj+=2){
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
              trans->beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;

       }
     } 

   swap1(image1, image, PIXEL **)

 }
  printf("\n");
}




void iterative_decoding_nonlinear(int level,int n_iter,double zoo)
{
  int rx,ry,rrx,rry,dx,dy;
  register int i,j,ii,jj,s;
  register PIXEL **imag,**imag1;
  double pixel;
  double z_factor;
  int width,height;
  static int first_time = 0;


  printf("\n");
  z_factor = zoo / (double) (1 << level);
  zooming(z_factor);
  width  = (int) rint(image_width  * z_factor / zoo); 
  height = (int) rint(image_height * z_factor / zoo); 

  if(first_time++ == 0)
     for(i=0;i< height;i++)
     for(j=0;j< width ;j++) 
        image[i][j] = 0;

  for(s=0; s < n_iter ; s++) {
     imag = image;
     imag1 = image1;

     if(level > 0)
         printf(" Decoding at low resolution (%dx%d) %d\r",width,height,s);
     else
         printf(" Iterative decoding ... %d\r",s);

     fflush(stdout);
     trans = &fractal_code;
     while(trans->next != NULL)  {
        trans = trans->next;
        rx=(int)rint(trans->rx);
        ry=(int)rint(trans->ry);
        dx=(int)rint(trans->dx);
        dy=(int)rint(trans->dy);
        rrx=(int)rint(trans->rrx);
        rry=(int)rint(trans->rry);

        switch(trans->sym_op) {     
         case IDENTITY   : 
            for(i=rx,ii=dx;i< rrx;i++,ii+=2)
            for(j=ry,jj=dy;j< rry;j++,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               imag1[i][j] = bound(0.5 +pixel * trans->alfa1 + pixel * pixel * trans->alfa2 +
                                                                      trans->beta);
               // if(imag1[i][j] > 127 ) imag1[i][j] = 127;
               // if(imag1[i][j] <  -128 ) imag1[i][j] = -128;
               // printf("image[i][j] = %d\n",image[i][j]);

               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case R_ROTATE90 :  
            for(j=rry-1,ii=dx;j>= ry;j--,ii+=2)
            for(i=rx,jj=dy;i< rrx;i++,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
                imag1[i][j] = bound(0.5 + pixel * trans->alfa1 + pixel * pixel * trans->alfa2 +
                                                                      trans->beta);
                // printf("image[i][j] = %d\n",image[i][j]);
                 // if(imag1[i][j] > 127 ) imag1[i][j] = 127;
                 // if(imag1[i][j] <  -128 ) imag1[i][j] = -128;
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }

            }
            break;
         case L_ROTATE90 :  
      for(j=ry,ii=dx;j< rry;j++,ii+=2) 
            for(i=rrx-1,jj=dy;i>= rx;i--,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
                imag1[i][j] = bound(0.5 +pixel * trans->alfa1 + pixel * pixel * trans->alfa2 +
                                                                      trans->beta);
                // printf("image[i][j] = %d\n",image[i][j]);
               //   if(imag1[i][j] > 127 ) imag1[i][j] = 127;
               // if(imag1[i][j] <  -128 ) imag1[i][j] = -128;
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }

            }
            break;
         case ROTATE180 :  
      for(i=rrx-1,ii=dx;i>= rx;i--,ii+=2) 
            for(j=rry-1,jj=dy;j>= ry;j--,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               imag1[i][j] = bound(0.5 +pixel * trans->alfa1 + pixel * pixel * trans->alfa2 +
                                                                      trans->beta);
               // printf("image[i][j] = %d\n",image[i][j]);
               //  if(imag1[i][j] > 127 ) imag1[i][j] = 127;
               // if(imag1[i][j] <  -128 ) imag1[i][j] = -128;
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case R_VERTICAL :  
      for(i=rx,ii=dx;i< rrx;i++,ii+=2) 
            for(j=rry-1,jj=dy;j>= ry;j--,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
                imag1[i][j] = bound(0.5 +pixel * trans->alfa1 + pixel * pixel * trans->alfa2 +
                                                                      trans->beta);
                // printf("image[i][j] = %d\n",image[i][j]);
               //   if(imag1[i][j] > 127 ) imag1[i][j] = 127;
               // if(imag1[i][j] <  -128 ) imag1[i][j] = -128;
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case R_HORIZONTAL: 
      for(i=rrx-1,ii=dx;i>= rx;i--,ii+=2)
            for(j=ry,jj=dy;j< rry;j++,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               imag1[i][j] = bound(0.5 +pixel * trans->alfa1 + pixel * pixel * trans->alfa2 +
                                                                      trans->beta);
               // printf("image[i][j] = %d\n",image[i][j]);
               //  if(imag1[i][j] > 127 ) imag1[i][j] = 127;
               // if(imag1[i][j] <  -128 ) imag1[i][j] = -128;
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case F_DIAGONAL :  
      for(j=ry,ii=dx;j< rry;j++,ii+=2)
            for(i=rx,jj=dy;i< rrx;i++,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               imag1[i][j] = bound(0.5 +pixel * trans->alfa1 + pixel * pixel * trans->alfa2 +
                                                                      trans->beta);
               // printf("image[i][j] = %d\n",image[i][j]);
               //  if(imag1[i][j] > 127 ) imag1[i][j] = 127;
               // if(imag1[i][j] <  -128 ) imag1[i][j] = -128;
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case S_DIAGONAL:   
      for(j=rry-1,ii=dx;j>= ry;j--,ii+=2) 
      for(i=rrx-1,jj=dy;i>= rx;i--,jj+=2){
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               imag1[i][j] = bound(0.5 +pixel * trans->alfa1 + pixel * pixel * trans->alfa2 +
                                                                      trans->beta);
               // printf("image[i][j] = %d\n",image[i][j]);
               //  if(imag1[i][j] > 127 ) imag1[i][j] = 127;
               // if(imag1[i][j] <  -128 ) imag1[i][j] = -128;
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;

       }
     } 

   swap1(image1, image, PIXEL **)

 }
  printf("\n");
}

void iterative_decoding_LumInv(int level,int n_iter,double zoo)
{
  int rx,ry,rrx,rry,dx,dy;
  register int i,j,ii,jj,s;
  register PIXEL **imag,**imag1;
  double pixel;
  double z_factor;
  int width,height;
  static int first_time = 0;


  printf("\n");
  z_factor = zoo / (double) (1 << level);
  zooming(z_factor);
  width  = (int) rint(image_width  * z_factor / zoo); 
  height = (int) rint(image_height * z_factor / zoo); 

  if(first_time++ == 0)
     for(i=0;i< height;i++)
     for(j=0;j< width ;j++) 
        image[i][j] = 128;

  for(s=0; s < n_iter ; s++) {
     imag = image;
     imag1 = image1;

     if(level > 0)
         printf(" Decoding at low resolution (%dx%d) %d\r",width,height,s);
     else
         printf(" Iterative decoding... %d\r",s);

     fflush(stdout);
     trans = &fractal_code;
     while(trans->next != NULL)  {
        trans = trans->next;
        rx=(int)rint(trans->rx);
        ry=(int)rint(trans->ry);
        dx=(int)rint(trans->dx);
        dy=(int)rint(trans->dy);
        rrx=(int)rint(trans->rrx);
        rry=(int)rint(trans->rry);

        switch(trans->sym_op) {     
         case IDENTITY   : 
            for(i=rx,ii=dx;i< rrx;i++,ii+=2)
            for(j=ry,jj=dy;j< rry;j++,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)+(imag[ii+1][jj] - pixel)+
                                    // (imag[ii][jj+1] - pixel)+(imag[ii+1][jj+1] - pixel))*trans->alfa + trans->rmean);
               // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)*trans->alfa + trans->rmean));
               imag1[i][j] = bound(0.5 + ((pixel)*trans->alfa + trans->rmean));
               // printf("%d\n",imag1[i][j]);

               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case R_ROTATE90 :  
            for(j=rry-1,ii=dx;j>= ry;j--,ii+=2)
            for(i=rx,jj=dy;i< rrx;i++,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)+(imag[ii+1][jj] - pixel)+
                                    // (imag[ii][jj+1] - pixel)+(imag[ii+1][jj+1] - pixel))*trans->alfa + trans->rmean);
               // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)*trans->alfa + trans->rmean));
               imag1[i][j] = bound(0.5 + ((pixel)*trans->alfa + trans->rmean));
               // printf("%d\n",imag1[i][j]);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }

            }
            break;
         case L_ROTATE90 :  
      for(j=ry,ii=dx;j< rry;j++,ii+=2) 
            for(i=rrx-1,jj=dy;i>= rx;i--,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)+(imag[ii+1][jj] - pixel)+
                                    // (imag[ii][jj+1] - pixel)+(imag[ii+1][jj+1] - pixel))*trans->alfa + trans->rmean);
               // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)*trans->alfa + trans->rmean));
               imag1[i][j] = bound(0.5 + ((pixel)*trans->alfa + trans->rmean));
               // printf("%d\n",imag1[i][j]);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }

            }
            break;
         case ROTATE180 :  
      for(i=rrx-1,ii=dx;i>= rx;i--,ii+=2) 
            for(j=rry-1,jj=dy;j>= ry;j--,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)+(imag[ii+1][jj] - pixel)+
                                    // (imag[ii][jj+1] - pixel)+(imag[ii+1][jj+1] - pixel))*trans->alfa + trans->rmean);
               // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)*trans->alfa + trans->rmean));
               imag1[i][j] = bound(0.5 + ((pixel)*trans->alfa + trans->rmean));
               // printf("%d\n",imag1[i][j]);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case R_VERTICAL :  
      for(i=rx,ii=dx;i< rrx;i++,ii+=2) 
            for(j=rry-1,jj=dy;j>= ry;j--,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)+(imag[ii+1][jj] - pixel)+
                                    // (imag[ii][jj+1] - pixel)+(imag[ii+1][jj+1] - pixel))*trans->alfa + trans->rmean);
               // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)*trans->alfa + trans->rmean));
               imag1[i][j] = bound(0.5 + ((pixel)*trans->alfa + trans->rmean));
               // printf("%d\n",imag1[i][j]);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case R_HORIZONTAL: 
      for(i=rrx-1,ii=dx;i>= rx;i--,ii+=2)
            for(j=ry,jj=dy;j< rry;j++,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)+(imag[ii+1][jj] - pixel)+
                                    // (imag[ii][jj+1] - pixel)+(imag[ii+1][jj+1] - pixel))*trans->alfa + trans->rmean);
               imag1[i][j] = bound(0.5 + ((pixel)*trans->alfa + trans->rmean));
               // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)*trans->alfa + trans->rmean));

               // printf("%d\n",imag1[i][j]);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case F_DIAGONAL :  
      for(j=ry,ii=dx;j< rry;j++,ii+=2)
            for(i=rx,jj=dy;i< rrx;i++,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)+(imag[ii+1][jj] - pixel)+
                                    // (imag[ii][jj+1] - pixel)+(imag[ii+1][jj+1] - pixel))*trans->alfa + trans->rmean);
               imag1[i][j] = bound(0.5 + ((pixel)*trans->alfa + trans->rmean));
               // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)*trans->alfa + trans->rmean));

               // printf("%d\n",imag1[i][j]);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case S_DIAGONAL:   
      for(j=rry-1,ii=dx;j>= ry;j--,ii+=2) 
      for(i=rrx-1,jj=dy;i>= rx;i--,jj+=2){
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)+(imag[ii+1][jj] - pixel)+
                                    // (imag[ii][jj+1] - pixel)+(imag[ii+1][jj+1] - pixel))*trans->alfa + trans->rmean);
               imag1[i][j] = bound(0.5 + ((pixel)*trans->alfa + trans->rmean));
               // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)*trans->alfa + trans->rmean));

               // printf("%d\n",imag1[i][j]);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;

       }
     } 

   swap1(image1, image, PIXEL **)

 }
   printf("\n");
}
void iterative_decoding_LumInv2(int level,int n_iter,double zoo)
{
  int rx,ry,rrx,rry,dx,dy;
  register int i,j,ii,jj,s;
  register PIXEL **imag,**imag1;
  double pixel;
  double z_factor;
  int width,height;
  static int first_time = 0;
  double beta  = 0.0;

  printf("\n");
  z_factor = zoo / (double) (1 << level);
  zooming(z_factor);
  width  = (int) rint(image_width  * z_factor / zoo); 
  height = (int) rint(image_height * z_factor / zoo); 

  if(first_time++ == 0)
     for(i=0;i< height;i++)
     for(j=0;j< width ;j++) 
        image[i][j] = 128;

  for(s=0; s < n_iter ; s++) {
     imag = image;
     imag1 = image1;

     if(level > 0)
          printf(" Decoding at low resolution (%dx%d) %d\r",width,height,s);
     else
          printf(" Iterative decoding... %d\r",s);

     fflush(stdout);
     trans = &fractal_code;
     while(trans->next != NULL)  {
        trans = trans->next;
        rx=(int)rint(trans->rx);
        ry=(int)rint(trans->ry);
        dx=(int)rint(trans->dx);
        dy=(int)rint(trans->dy);
        rrx=(int)rint(trans->rrx);
        rry=(int)rint(trans->rry);

        switch(trans->sym_op) {     
         case IDENTITY   : 
            for(i=rx,ii=dx;i< rrx;i++,ii+=2)
            for(j=ry,jj=dy;j< rry;j++,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + beta);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case R_ROTATE90 :  
            for(j=rry-1,ii=dx;j>= ry;j--,ii+=2)
            for(i=rx,jj=dy;i< rrx;i++,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + beta);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }

            }
            break;
         case L_ROTATE90 :  
      for(j=ry,ii=dx;j< rry;j++,ii+=2) 
            for(i=rrx-1,jj=dy;i>= rx;i--,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + beta);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }

            }
            break;
         case ROTATE180 :  
      for(i=rrx-1,ii=dx;i>= rx;i--,ii+=2) 
            for(j=rry-1,jj=dy;j>= ry;j--,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + beta);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case R_VERTICAL :  
      for(i=rx,ii=dx;i< rrx;i++,ii+=2) 
            for(j=rry-1,jj=dy;j>= ry;j--,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + beta);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case R_HORIZONTAL: 
      for(i=rrx-1,ii=dx;i>= rx;i--,ii+=2)
            for(j=ry,jj=dy;j< rry;j++,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + beta);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case F_DIAGONAL :  
      for(j=ry,ii=dx;j< rry;j++,ii+=2)
            for(i=rx,jj=dy;i< rrx;i++,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + beta);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case S_DIAGONAL:   
      for(j=rry-1,ii=dx;j>= ry;j--,ii+=2) 
      for(i=rrx-1,jj=dy;i>= rx;i--,jj+=2){
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + beta);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;

       }
     } 

   swap1(image1, image, PIXEL **)

 }
  printf("\n");
}

void iterative_decoding_new_init_image(int level,int n_iter,double zoo)
{
  int rx,ry,rrx,rry,dx,dy;
  register int i,j,ii,jj,s;
  register PIXEL **imag,**imag1;
  double pixel;
  double z_factor;
  int width,height;
  static int first_time = 0;


  printf("\n");
  z_factor = zoo / (double) (1 << level);
  zooming(z_factor);
  width  = (int) rint(image_width  * z_factor / zoo); 
  height = (int) rint(image_height * z_factor / zoo); 
  // Mat startImg;
  // startImg = imread("./filter/noise_free.png",CV_LOAD_IMAGE_GRAYSCALE);
  // if(first_time++ == 0)
     // for(i=0;i< height;i++)
     // for(j=0;j< width ;j++){ 
     //     printf("image[%d][%d] = %u\n",i,j,(unsigned int)image[i][j]);
     //    // image[i][j] = startImg.at<uchar>(i,j);  
     //    // image[i][j] = 128;
     //  }


  for(s=0; s < n_iter ; s++) {
     imag = image;
     imag1 = image1;

     if(level > 0)
         printf(" Decoding at low resolution (%dx%d) %d\r",width,height,s);
     else
         printf(" Iterative decoding... %d\r",s);

     fflush(stdout);
     trans = &fractal_code;
     while(trans->next != NULL)  {
        trans = trans->next;
        rx=(int)rint(trans->rx);
        ry=(int)rint(trans->ry);
        dx=(int)rint(trans->dx);
        dy=(int)rint(trans->dy);
        rrx=(int)rint(trans->rrx);
        rry=(int)rint(trans->rry);

        switch(trans->sym_op) {     
         case IDENTITY   : 
            for(i=rx,ii=dx;i< rrx;i++,ii+=2)
            for(j=ry,jj=dy;j< rry;j++,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
               // printf("%d\n",imag1[i][j]);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case R_ROTATE90 :  
            for(j=rry-1,ii=dx;j>= ry;j--,ii+=2)
            for(i=rx,jj=dy;i< rrx;i++,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
               // printf("%d\n",imag1[i][j]);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }

            }
            break;
         case L_ROTATE90 :  
      for(j=ry,ii=dx;j< rry;j++,ii+=2) 
            for(i=rrx-1,jj=dy;i>= rx;i--,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
               // printf("%d\n",imag1[i][j]);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }

            }
            break;
         case ROTATE180 :  
      for(i=rrx-1,ii=dx;i>= rx;i--,ii+=2) 
            for(j=rry-1,jj=dy;j>= ry;j--,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
               // printf("%d\n",imag1[i][j]);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case R_VERTICAL :  
      for(i=rx,ii=dx;i< rrx;i++,ii+=2) 
            for(j=rry-1,jj=dy;j>= ry;j--,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
               // printf("%d\n",imag1[i][j]);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case R_HORIZONTAL: 
      for(i=rrx-1,ii=dx;i>= rx;i--,ii+=2)
            for(j=ry,jj=dy;j< rry;j++,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
               // printf("%d\n",imag1[i][j]);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case F_DIAGONAL :  
      for(j=ry,ii=dx;j< rry;j++,ii+=2)
            for(i=rx,jj=dy;i< rrx;i++,jj+=2) {
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
               // printf("%d\n",imag1[i][j]);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;
         case S_DIAGONAL:   
      for(j=rry-1,ii=dx;j>= ry;j--,ii+=2) 
      for(i=rrx-1,jj=dy;i>= rx;i--,jj+=2){
               pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
               // printf("%d\n",imag1[i][j]);
               if (isColor) {
                 image_uch[i][j] = trans->um;
                 image_vch[i][j] = trans->vm;
               }
            }
            break;

       }
     } 

   swap1(image1, image, PIXEL **)

 }
  printf("\n");
}


void iterative_decoding(int level,int n_iter,double zoo)
{
  int rx,ry,rrx,rry,dx,dy;
  register int i,j,ii,jj,s;
  register PIXEL **imag,**imag1;
  double pixel;
  double z_factor;
  int width,height;
  static int first_time = 0;


  printf("\n");
  z_factor = zoo / (double) (1 << level);
  zooming(z_factor);
  width  = (int) rint(image_width  * z_factor / zoo); 
  height = (int) rint(image_height * z_factor / zoo); 
  // Mat startImg;
  // startImg = imread("./filter/noise_free.png",CV_LOAD_IMAGE_GRAYSCALE);
  if(first_time++ == 0)
     for(i=0;i< height;i++)
     for(j=0;j< width ;j++){ 
         // printf("image[%d][%d] = %u\n",i,j,image[i][j]);
        // image[i][j] = startImg.at<uchar>(i,j);  
        image[i][j] = 128;
      }


  for(s=0; s < n_iter ; s++) {
   imag = image;
   imag1 = image1;

   if(level > 0)
       printf(" Decoding at low resolution (%dx%d) %d\r",width,height,s);
   else
       printf(" Iterative decoding... %d\r",s);

   fflush(stdout);
   trans = &fractal_code;
   while(trans->next != NULL)  {
      trans = trans->next;
      rx=(int)rint(trans->rx);
      ry=(int)rint(trans->ry);
      dx=(int)rint(trans->dx);
      dy=(int)rint(trans->dy);
      rrx=(int)rint(trans->rrx);
      rry=(int)rint(trans->rry);

      switch(trans->sym_op) {     
       case IDENTITY   : 
          for(i=rx,ii=dx;i< rrx;i++,ii+=2)
          for(j=ry,jj=dy;j< rry;j++,jj+=2) {
             pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                       imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
             imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
             // printf("%d\n",imag1[i][j]);
             if (isColor) {
               image_uch[i][j] = trans->um;
               image_vch[i][j] = trans->vm;
             }
          }
          break;
       case R_ROTATE90 :  
          for(j=rry-1,ii=dx;j>= ry;j--,ii+=2)
          for(i=rx,jj=dy;i< rrx;i++,jj+=2) {
             pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                       imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
             imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
             // printf("%d\n",imag1[i][j]);
             if (isColor) {
               image_uch[i][j] = trans->um;
               image_vch[i][j] = trans->vm;
             }

          }
          break;
       case L_ROTATE90 :  
    for(j=ry,ii=dx;j< rry;j++,ii+=2) 
          for(i=rrx-1,jj=dy;i>= rx;i--,jj+=2) {
             pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                       imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
             imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
             // printf("%d\n",imag1[i][j]);
             if (isColor) {
               image_uch[i][j] = trans->um;
               image_vch[i][j] = trans->vm;
             }

          }
          break;
       case ROTATE180 :  
    for(i=rrx-1,ii=dx;i>= rx;i--,ii+=2) 
          for(j=rry-1,jj=dy;j>= ry;j--,jj+=2) {
             pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                       imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
             imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
             // printf("%d\n",imag1[i][j]);
             if (isColor) {
               image_uch[i][j] = trans->um;
               image_vch[i][j] = trans->vm;
             }
          }
          break;
       case R_VERTICAL :  
    for(i=rx,ii=dx;i< rrx;i++,ii+=2) 
          for(j=rry-1,jj=dy;j>= ry;j--,jj+=2) {
             pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                       imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
             imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
             // printf("%d\n",imag1[i][j]);
             if (isColor) {
               image_uch[i][j] = trans->um;
               image_vch[i][j] = trans->vm;
             }
          }
          break;
       case R_HORIZONTAL: 
    for(i=rrx-1,ii=dx;i>= rx;i--,ii+=2)
          for(j=ry,jj=dy;j< rry;j++,jj+=2) {
             pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                       imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
             imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
             // printf("%d\n",imag1[i][j]);
             if (isColor) {
               image_uch[i][j] = trans->um;
               image_vch[i][j] = trans->vm;
             }
          }
          break;
       case F_DIAGONAL :  
    for(j=ry,ii=dx;j< rry;j++,ii+=2)
          for(i=rx,jj=dy;i< rrx;i++,jj+=2) {
             pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                       imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
             imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
             // printf("%d\n",imag1[i][j]);
             if (isColor) {
               image_uch[i][j] = trans->um;
               image_vch[i][j] = trans->vm;
             }
          }
          break;
       case S_DIAGONAL:   
    for(j=rry-1,ii=dx;j>= ry;j--,ii+=2) 
    for(i=rrx-1,jj=dy;i>= rx;i--,jj+=2){
             pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                       imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
             imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
             // printf("%d\n",imag1[i][j]);
             if (isColor) {
               image_uch[i][j] = trans->um;
               image_vch[i][j] = trans->vm;
             }
          }
          break;

     }
   } 
   swap1(image1, image, PIXEL **)
  }
  printf("\n");
}

void piramidal_decoding_LumInv(int level)
{
  int rx,ry,rrx,rry,dx,dy;
  register int i,j,ii,jj;
  register PIXEL **imag,**imag1;
  double pixel;

  if(level < 1)
     return;

  zooming(2.0);   /* Increase resolution */

  imag = image;
  imag1 = image1;

  printf(" Increasing resolution... \r");
  fflush(stdout);
  trans = &fractal_code;

  while(trans->next != NULL)  {
     trans = trans->next;

     rx=(int)rint(trans->rx);
     ry=(int)rint(trans->ry);
     dx=(int)rint(trans->dx);
     dy=(int)rint(trans->dy);
     rrx=(int)rint(trans->rrx);
     rry=(int)rint(trans->rry);

     switch(trans->sym_op) {
        case IDENTITY   :  
     for(i=rx,ii=dx >> 1;i< rrx;i++,ii++)
           for(j=ry,jj=dy >> 1;j< rry;j++,jj++) {
              pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)+(imag[ii+1][jj] - pixel)+
                                    // (imag[ii][jj+1] - pixel)+(imag[ii+1][jj+1] - pixel))*trans->alfa + trans->rmean);
              imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)*trans->alfa + trans->rmean));
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case R_ROTATE90 :  
     for(j=rry-1,ii=dx >> 1;j>= ry;j--,ii++)
           for(i=rx,jj=dy >> 1;i< rrx;i++,jj++) {
             pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)+(imag[ii+1][jj] - pixel)+
                                    // (imag[ii][jj+1] - pixel)+(imag[ii+1][jj+1] - pixel))*trans->alfa + trans->rmean);
             // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)*trans->alfa + trans->rmean));
             imag1[i][j] = bound(0.5 + ((pixel)*trans->alfa + trans->rmean));
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case L_ROTATE90 :  
     for(j=ry,ii=dx >> 1;j< rry;j++,ii++)
           for(i=rrx-1,jj=dy >> 1;i>= rx;i--,jj++) {
              pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)+(imag[ii+1][jj] - pixel)+
                                    // (imag[ii][jj+1] - pixel)+(imag[ii+1][jj+1] - pixel))*trans->alfa + trans->rmean);
              // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)*trans->alfa + trans->rmean));
              imag1[i][j] = bound(0.5 + ((pixel)*trans->alfa + trans->rmean));
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case ROTATE180 :  
     for(i=rrx-1,ii=dx >> 1;i>= rx;i--,ii++)
           for(j=rry-1,jj=dy >> 1;j>= ry;j--,jj++) {
            pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)+(imag[ii+1][jj] - pixel)+
                                    // (imag[ii][jj+1] - pixel)+(imag[ii+1][jj+1] - pixel))*trans->alfa + trans->rmean);
            // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)*trans->alfa + trans->rmean));
            imag1[i][j] = bound(0.5 + ((pixel)*trans->alfa + trans->rmean));
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case R_VERTICAL :  
     for(i=rx,ii=dx >> 1;i< rrx;i++,ii++)
           for(j=rry-1,jj=dy >> 1;j>= ry;j--,jj++) {
              pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)+(imag[ii+1][jj] - pixel)+
                                    // (imag[ii][jj+1] - pixel)+(imag[ii+1][jj+1] - pixel))*trans->alfa + trans->rmean);
              // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)*trans->alfa + trans->rmean));
              imag1[i][j] = bound(0.5 + ((pixel)*trans->alfa + trans->rmean));
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case R_HORIZONTAL: 
     for(i=rrx-1,ii=dx >> 1;i>= rx;i--,ii++)
           for(j=ry,jj=dy >> 1;j< rry;j++,jj++) {
             pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)+(imag[ii+1][jj] - pixel)+
                                    // (imag[ii][jj+1] - pixel)+(imag[ii+1][jj+1] - pixel))*trans->alfa + trans->rmean);
             // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)*trans->alfa + trans->rmean));
             imag1[i][j] = bound(0.5 + ((pixel)*trans->alfa + trans->rmean));
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case F_DIAGONAL :  
     for(j=ry,ii=dx >> 1;j< rry;j++,ii++)
           for(i=rx,jj=dy >> 1;i< rrx;i++,jj++) {
            pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)+(imag[ii+1][jj] - pixel)+
                                    // (imag[ii][jj+1] - pixel)+(imag[ii+1][jj+1] - pixel))*trans->alfa + trans->rmean);
            // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)*trans->alfa + trans->rmean));
            imag1[i][j] = bound(0.5 + ((pixel)*trans->alfa + trans->rmean));
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case S_DIAGONAL:   
     for(j=rry-1,ii=dx >> 1;j>= ry;j--,ii++)
           for(i=rrx-1,jj=dy >> 1;i>= rx;i--,jj++){
             pixel = (double)(imag[ii][jj]+imag[ii+1][jj] +
                                         imag[ii][jj+1]+imag[ii+1][jj+1])/4.0;
               // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)+(imag[ii+1][jj] - pixel)+
                                    // (imag[ii][jj+1] - pixel)+(imag[ii+1][jj+1] - pixel))*trans->alfa + trans->rmean);
             // imag1[i][j] = bound(0.5 + ((imag[ii][jj] - pixel)*trans->alfa + trans->rmean));
             imag1[i][j] = bound(0.5 + ((pixel)*trans->alfa + trans->rmean));
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
     }

  } 
  swap1(image1, image, PIXEL **)

  if(level > 1) {
    if(quality)
       iterative_decoding_LumInv(0,2,1.0); 
    piramidal_decoding_LumInv(level-1);
  }

}

void piramidal_decoding_testing(int level)
{
  int rx,ry,rrx,rry,dx,dy;
  register int i,j,ii,jj;
  register PIXEL **imag,**imag1;
  double pixel;

  if(level < 1)
     return;

  zooming(2.0);   /* Increase resolution */

  imag = image;
  imag1 = image1;

  printf(" Increasing resolution... \r");
  fflush(stdout);
  trans = &fractal_code;

  while(trans->next != NULL)  {
     trans = trans->next;

     rx=(int)rint(trans->rx);
     ry=(int)rint(trans->ry);
     dx=(int)rint(trans->dx);
     dy=(int)rint(trans->dy);
     rrx=(int)rint(trans->rrx);
     rry=(int)rint(trans->rry);

     switch(trans->sym_op) {
        case IDENTITY   :  
     for(i=rx,ii=dx >> 1;i< rrx;i++,ii++)
           for(j=ry,jj=dy >> 1;j< rry;j++,jj++) {
              pixel = (double)imag[ii][jj];
              trans->beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case R_ROTATE90 :  
     for(j=rry-1,ii=dx >> 1;j>= ry;j--,ii++)
           for(i=rx,jj=dy >> 1;i< rrx;i++,jj++) {
              pixel = (double)imag[ii][jj];
              trans->beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case L_ROTATE90 :  
     for(j=ry,ii=dx >> 1;j< rry;j++,ii++)
           for(i=rrx-1,jj=dy >> 1;i>= rx;i--,jj++) {
              pixel = (double)imag[ii][jj];

              trans->beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case ROTATE180 :  
     for(i=rrx-1,ii=dx >> 1;i>= rx;i--,ii++)
           for(j=rry-1,jj=dy >> 1;j>= ry;j--,jj++) {
              pixel = (double)imag[ii][jj];
              trans->beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case R_VERTICAL :  
     for(i=rx,ii=dx >> 1;i< rrx;i++,ii++)
           for(j=rry-1,jj=dy >> 1;j>= ry;j--,jj++) {
              pixel = (double)imag[ii][jj];
              trans->beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case R_HORIZONTAL: 
     for(i=rrx-1,ii=dx >> 1;i>= rx;i--,ii++)
           for(j=ry,jj=dy >> 1;j< rry;j++,jj++) {
              pixel = (double)imag[ii][jj];
              trans->beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case F_DIAGONAL :  
     for(j=ry,ii=dx >> 1;j< rry;j++,ii++)
           for(i=rx,jj=dy >> 1;i< rrx;i++,jj++) {
              pixel = (double)imag[ii][jj];
              
              trans->beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case S_DIAGONAL:   
     for(j=rry-1,ii=dx >> 1;j>= ry;j--,ii++)
           for(i=rrx-1,jj=dy >> 1;i>= rx;i--,jj++){
              pixel = (double)imag[ii][jj];
              trans->beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
     }

  } 
  swap1(image1, image, PIXEL **)

  if(level > 1) {
    if(quality)
       iterative_decoding_testing(0,2,1.0); 
    piramidal_decoding_testing(level-1);
  }

}

void piramidal_decoding_nonlinear(int level)
{
  int rx,ry,rrx,rry,dx,dy;
  register int i,j,ii,jj;
  register PIXEL **imag,**imag1;
  double pixel;

  if(level < 1)
     return;

  zooming(2.0);   /* Increase resolution */

  imag = image;
  imag1 = image1;

  printf(" Increasing resolution... \r");
  fflush(stdout);
  trans = &fractal_code;

  while(trans->next != NULL)  {
     trans = trans->next;

     rx=(int)rint(trans->rx);
     ry=(int)rint(trans->ry);
     dx=(int)rint(trans->dx);
     dy=(int)rint(trans->dy);
     rrx=(int)rint(trans->rrx);
     rry=(int)rint(trans->rry);

     switch(trans->sym_op) {
        case IDENTITY   :  
     for(i=rx,ii=dx >> 1;i< rrx;i++,ii++)
           for(j=ry,jj=dy >> 1;j< rry;j++,jj++) {
              pixel = (double)imag[ii][jj];
               imag1[i][j] = bound(0.5 +pixel * trans->alfa1 + pixel * pixel *trans->alfa2 +
                                                                      trans->beta);
               //  if(imag1[i][j] > 127 ) imag1[i][j] = 127;
               // if(imag1[i][j] <  -128 ) imag1[i][j] = -128;
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case R_ROTATE90 :  
     for(j=rry-1,ii=dx >> 1;j>= ry;j--,ii++)
           for(i=rx,jj=dy >> 1;i< rrx;i++,jj++) {
              pixel = (double)imag[ii][jj];
               imag1[i][j] = bound(0.5 +pixel * trans->alfa1 + pixel * pixel *trans->alfa2 +
                                                                      trans->beta);
               //  if(imag1[i][j] > 127 ) imag1[i][j] = 127;
               // if(imag1[i][j] <  -128 ) imag1[i][j] = -128;
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case L_ROTATE90 :  
     for(j=ry,ii=dx >> 1;j< rry;j++,ii++)
           for(i=rrx-1,jj=dy >> 1;i>= rx;i--,jj++) {
              pixel = (double)imag[ii][jj];
               imag1[i][j] = bound(0.5 + pixel * trans->alfa1 + pixel * pixel *trans->alfa2 +
                                                                      trans->beta);
               //  if(imag1[i][j] > 127 ) imag1[i][j] = 127;
               // if(imag1[i][j] <  -128 ) imag1[i][j] = -128;
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case ROTATE180 :  
     for(i=rrx-1,ii=dx >> 1;i>= rx;i--,ii++)
           for(j=rry-1,jj=dy >> 1;j>= ry;j--,jj++) {
              pixel = (double)imag[ii][jj];
              imag1[i][j] = bound(0.5 +pixel * trans->alfa1 + pixel * pixel *trans->alfa2 +
                                                                      trans->beta);
               // if(imag1[i][j] > 127 ) imag1[i][j] = 127;
               // if(imag1[i][j] <  -128 ) imag1[i][j] = -128;
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case R_VERTICAL :  
     for(i=rx,ii=dx >> 1;i< rrx;i++,ii++)
           for(j=rry-1,jj=dy >> 1;j>= ry;j--,jj++) {
              pixel = (double)imag[ii][jj];
               imag1[i][j] = bound(0.5 +pixel * trans->alfa1 + pixel * pixel *trans->alfa2 +
                                                                      trans->beta);
               //  if(imag1[i][j] > 127 ) imag1[i][j] = 127;
               // if(imag1[i][j] <  -128 ) imag1[i][j] = -128;
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case R_HORIZONTAL: 
     for(i=rrx-1,ii=dx >> 1;i>= rx;i--,ii++)
           for(j=ry,jj=dy >> 1;j< rry;j++,jj++) {
              pixel = (double)imag[ii][jj];
               imag1[i][j] = bound(0.5 +pixel * trans->alfa1 + pixel * pixel *trans->alfa2 +
                                                                      trans->beta);
               //  if(imag1[i][j] > 127 ) imag1[i][j] = 127;
               // if(imag1[i][j] <  -128 ) imag1[i][j] = -128;
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case F_DIAGONAL :  
     for(j=ry,ii=dx >> 1;j< rry;j++,ii++)
           for(i=rx,jj=dy >> 1;i< rrx;i++,jj++) {
              pixel = (double)imag[ii][jj];
               imag1[i][j] = bound(0.5 +pixel * trans->alfa1 + pixel * pixel *trans->alfa2 +
                                                                      trans->beta);
               //  if(imag1[i][j] > 127 ) imag1[i][j] = 127;
               // if(imag1[i][j] <  -128 ) imag1[i][j] = -128;
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case S_DIAGONAL:   
     for(j=rry-1,ii=dx >> 1;j>= ry;j--,ii++)
           for(i=rrx-1,jj=dy >> 1;i>= rx;i--,jj++){
              pixel = (double)imag[ii][jj];
             imag1[i][j] = bound(0.5 + pixel * trans->alfa1 + pixel * pixel *trans->alfa2 +
                                                                      trans->beta);
              // if(imag1[i][j] > 127 ) imag1[i][j] = 127;
              //  if(imag1[i][j] <  -128 ) imag1[i][j] = -128;
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
     }

  } 
  swap1(image1, image, PIXEL **)

  if(level > 1) {
    if(quality)
       iterative_decoding_nonlinear(0,2,1.0); 
    piramidal_decoding_nonlinear(level-1);
  }

}
void piramidal_decoding_LumInv2(int level)
{
  int rx,ry,rrx,rry,dx,dy;
  register int i,j,ii,jj;
  register PIXEL **imag,**imag1;
  double pixel;
  double beta;
  if(level < 1)
     return;

  zooming(2.0);   /* Increase resolution */

  imag = image;
  imag1 = image1;

  printf(" Increasing resolution... \r");
  fflush(stdout);
  trans = &fractal_code;

  while(trans->next != NULL)  {
     trans = trans->next;

     rx=(int)rint(trans->rx);
     ry=(int)rint(trans->ry);
     dx=(int)rint(trans->dx);
     dy=(int)rint(trans->dy);
     rrx=(int)rint(trans->rrx);
     rry=(int)rint(trans->rry);

     switch(trans->sym_op) {
        case IDENTITY   :  
     for(i=rx,ii=dx >> 1;i< rrx;i++,ii++)
           for(j=ry,jj=dy >> 1;j< rry;j++,jj++) {
              pixel = (double)imag[ii][jj];
              beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + beta);
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case R_ROTATE90 :  
     for(j=rry-1,ii=dx >> 1;j>= ry;j--,ii++)
           for(i=rx,jj=dy >> 1;i< rrx;i++,jj++) {
              pixel = (double)imag[ii][jj];
              beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + beta);
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case L_ROTATE90 :  
     for(j=ry,ii=dx >> 1;j< rry;j++,ii++)
           for(i=rrx-1,jj=dy >> 1;i>= rx;i--,jj++) {
              pixel = (double)imag[ii][jj];
              beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + beta);
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case ROTATE180 :  
     for(i=rrx-1,ii=dx >> 1;i>= rx;i--,ii++)
           for(j=rry-1,jj=dy >> 1;j>= ry;j--,jj++) {
              pixel = (double)imag[ii][jj];
              beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + beta);
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case R_VERTICAL :  
     for(i=rx,ii=dx >> 1;i< rrx;i++,ii++)
           for(j=rry-1,jj=dy >> 1;j>= ry;j--,jj++) {
              pixel = (double)imag[ii][jj];
              beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + beta);
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case R_HORIZONTAL: 
     for(i=rrx-1,ii=dx >> 1;i>= rx;i--,ii++)
           for(j=ry,jj=dy >> 1;j< rry;j++,jj++) {
              pixel = (double)imag[ii][jj];
              beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + beta);
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case F_DIAGONAL :  
     for(j=ry,ii=dx >> 1;j< rry;j++,ii++)
           for(i=rx,jj=dy >> 1;i< rrx;i++,jj++) {
              pixel = (double)imag[ii][jj];
              beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + beta);
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case S_DIAGONAL:   
     for(j=rry-1,ii=dx >> 1;j>= ry;j--,ii++)
           for(i=rrx-1,jj=dy >> 1;i>= rx;i--,jj++){
              pixel = (double)imag[ii][jj];
              beta = trans->rmean - trans->alfa * pixel;
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + beta);
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
     }

  } 
  swap1(image1, image, PIXEL **)

  if(level > 1) {
    if(quality)
       iterative_decoding(0,2,1.0); 
    piramidal_decoding(level-1);
  }

}
void piramidal_decoding(int level)
{
  printf("using piramidal decoding\n");
  int rx,ry,rrx,rry,dx,dy;
  register int i,j,ii,jj;
  register PIXEL **imag,**imag1;
  double pixel;

  if(level < 1)
     return;

  zooming(2.0);   /* Increase resolution */

  imag = image;
  imag1 = image1;

  printf(" Increasing resolution... \r");
  fflush(stdout);
  trans = &fractal_code;

  while(trans->next != NULL)  {
     trans = trans->next;

     rx=(int)rint(trans->rx);
     ry=(int)rint(trans->ry);
     dx=(int)rint(trans->dx);
     dy=(int)rint(trans->dy);
     rrx=(int)rint(trans->rrx);
     rry=(int)rint(trans->rry);

     switch(trans->sym_op) {
        case IDENTITY   :  
	   for(i=rx,ii=dx >> 1;i< rrx;i++,ii++)
           for(j=ry,jj=dy >> 1;j< rry;j++,jj++) {
              pixel = (double)imag[ii][jj];
              imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case R_ROTATE90 :  
	   for(j=rry-1,ii=dx >> 1;j>= ry;j--,ii++)
           for(i=rx,jj=dy >> 1;i< rrx;i++,jj++) {
              pixel = (double)imag[ii][jj];
              imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case L_ROTATE90 :  
	   for(j=ry,ii=dx >> 1;j< rry;j++,ii++)
           for(i=rrx-1,jj=dy >> 1;i>= rx;i--,jj++) {
              pixel = (double)imag[ii][jj];
              imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case ROTATE180 :  
	   for(i=rrx-1,ii=dx >> 1;i>= rx;i--,ii++)
           for(j=rry-1,jj=dy >> 1;j>= ry;j--,jj++) {
              pixel = (double)imag[ii][jj];
              imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case R_VERTICAL :  
	   for(i=rx,ii=dx >> 1;i< rrx;i++,ii++)
           for(j=rry-1,jj=dy >> 1;j>= ry;j--,jj++) {
              pixel = (double)imag[ii][jj];
              imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case R_HORIZONTAL: 
	   for(i=rrx-1,ii=dx >> 1;i>= rx;i--,ii++)
           for(j=ry,jj=dy >> 1;j< rry;j++,jj++) {
              pixel = (double)imag[ii][jj];
              imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case F_DIAGONAL :  
	   for(j=ry,ii=dx >> 1;j< rry;j++,ii++)
           for(i=rx,jj=dy >> 1;i< rrx;i++,jj++) {
              pixel = (double)imag[ii][jj];
              imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
        case S_DIAGONAL:   
	   for(j=rry-1,ii=dx >> 1;j>= ry;j--,ii++)
           for(i=rrx-1,jj=dy >> 1;i>= rx;i--,jj++){
              pixel = (double)imag[ii][jj];
              imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
              if (isColor) {
                image_uch[i][j] = trans->um;
                image_vch[i][j] = trans->vm;
              }
           }
           break;
     }

  } 
  swap1(image1, image, PIXEL **)

  if(level > 1) {
    if(quality)
       iterative_decoding(0,2,1.0); 
    piramidal_decoding(level-1);
  }

}


void smooth_image()
{
    double pixel1, pixel2;
    int i,j;
    int w1,w2;
    int rx,ry,rrx,rry;

    printf("\n Postprocessing ...");
    fflush(stdout);

    trans = &fractal_code;
    while (trans->next != NULL) {
       trans = trans->next;

       rx=(int)rint(trans->rx);
       ry=(int)rint(trans->ry);
       rrx=(int)rint(trans->rrx);
       rry=(int)rint(trans->rry);

       if (rx == 0 || ry == 0 || (int)trans->size == 1)
           continue;

       if (trans->size == min_size) { 
          w1 = 5;
          w2 = 1;
       } 
       else {
          w1 = 2;
          w2 = 1;
       }

       for (i=rx; i<rrx; ++i) {
          pixel1 = image[i][ry];
          pixel2 = image[i][ry-1];
          image[i][ry] = (w1*pixel1 + w2*pixel2)/(w1+w2);
          image[i][ry-1] = (w2*pixel1 + w1*pixel2)/(w1+w2);
       }

       for (j=ry; j<rry; ++j) {
          pixel1 = image[rx][j];
          pixel2 = image[rx-1][j];
          image[rx][j] = (w1*pixel1 + w2*pixel2)/(w1+w2);
          image[rx-1][j] = (w2*pixel1 + w1*pixel2)/(w1+w2);
       }
   }
   printf("done \n");
   fflush(stdout);
}



