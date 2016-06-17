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
#include <vector>

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>



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
  isNonlinear    = (int)unpack(1,input);
  printf("Nonlinear : %s\n",isNonlinear ? "True" :"False");

  if(isNonlinear){
  N_BITALFA1 =  (int)unpack(4,input);
  N_BITALFA2 =  (int)unpack(4,input);
  N_BITBETA2 =  (int)unpack(4,input);
  
  }else{
   N_BITALFA = (int)unpack(4,input);
   N_BITBETA = (int)unpack(4,input);
  }
  min_size       = (int)unpack(7,input);
  max_size       = (int)unpack(7,input);
  SHIFT          = (int)unpack(6,input);
  image_width    = (int)unpack(12,input);
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
  if(isNonlinear){
    read_transformations_nonlinear(0,0,virtual_size);
  }else
    read_transformations(0,0,virtual_size);
  printf("Tranforms = %d\n",transforms);
  printf("done\n");
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
     }else
      iterative_decoding(lev,iterations,zoom);
     if(isNonlinear){ 
       piramidal_decoding_nonlinear(lev); /* Increase resolution      */
     }else
       piramidal_decoding(lev);                     
     if(quality){
        if(isNonlinear)
          iterative_decoding_nonlinear(0,2,1.0); 
        else 
          iterative_decoding(0,2,1.0);     
     }    
  } else{
      if(isNonlinear)
        iterative_decoding_nonlinear(0,iterations,zoom);
      else
        iterative_decoding(0,iterations,zoom); 
    }
  

  if(postproc) 
     smooth_image(); 

  if( isColor ){ // Conver to RGB
    std::vector<cv::Mat> yuvChannels;
    cv::Mat ych(image_width,image_height, CV_8U);
    cv::Mat uch(image_width,image_height, CV_8U);
    cv::Mat vch(image_width,image_height, CV_8U);
    yuvChannels.push_back(ych);
    yuvChannels.push_back(uch);
    yuvChannels.push_back(vch);

    for (int iii = 0; iii < image_width; iii++) {
      for (int jjj = 0; jjj < image_height; jjj++) {
        yuvChannels[0].at<uchar>(jjj, iii) = image[jjj][iii];
        yuvChannels[1].at<uchar>(jjj, iii) = image_uch[jjj][iii];
        yuvChannels[2].at<uchar>(jjj, iii) = image_vch[jjj][iii];
      }
    }
    cv::merge(yuvChannels, ych);
    cv::cvtColor(ych, ych, CV_YUV2BGR);
    printf("Writing... output.dec.ppm");
    cv::imwrite("output.dec.ppm", ych);
    
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
    if(raw_format)  
       writeimage_raw(fileout, image,image_width,image_height);
    else
      writeimage_pgm(fileout, image,image_width,image_height);
     
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
       alfa1 = (qalfa1 - 15) * 0.0001;
       alfa2 = (qalfa2 - 15) * 0.0001;
      

      
      /* Compute beta from the quantized value */
      beta = (qbeta - 128);
      // alfa2 = (double) qalfa2 / (double)(1 << N_BITALFA2) * ( MAX_ALFA) ;
      // beta = (double)qbeta/(double)((1 << N_BITBETA2)-1)*((1.0+fabs(alfa1))*255);
      // if (alfa1 > 0.0) beta  -= alfa1 * 255;
         
      trans->alfa1 = alfa1;
      trans->alfa2 = alfa2;
      trans->beta = beta;
      if(isColor){
        trans->um = um;
        trans->vm = vm;
      }
     if(qalfa1 != zeroalfa) {      
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
        //  printf("%d %d\n", ddx, ddy);
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
               imag1[i][j] = bound(0.5 + pixel * trans->alfa1 + pixel * pixel * trans->alfa2 +
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
                imag1[i][j] = bound(0.5 + pixel * trans->alfa1 + pixel * pixel * trans->alfa2 +
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
               imag1[i][j] = bound(0.5 + pixel * trans->alfa1 + pixel * pixel * trans->alfa2 +
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
                imag1[i][j] = bound(0.5 + pixel * trans->alfa1 + pixel * pixel * trans->alfa2 +
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
               imag1[i][j] = bound(0.5 + pixel * trans->alfa1 + pixel * pixel * trans->alfa2 +
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
               imag1[i][j] = bound(0.5 + pixel * trans->alfa1 + pixel * pixel * trans->alfa2 +
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
               imag1[i][j] = bound(0.5 + pixel * trans->alfa1 + pixel * pixel * trans->alfa2 +
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
               imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
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
        case R_ROTATE90 :  
     for(j=rry-1,ii=dx >> 1;j>= ry;j--,ii++)
           for(i=rx,jj=dy >> 1;i< rrx;i++,jj++) {
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
              imag1[i][j] = bound(0.5 + pixel * trans->alfa1 + pixel * pixel *trans->alfa2 +
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
        case R_HORIZONTAL: 
     for(i=rrx-1,ii=dx >> 1;i>= rx;i--,ii++)
           for(j=ry,jj=dy >> 1;j< rry;j++,jj++) {
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
        case F_DIAGONAL :  
     for(j=ry,ii=dx >> 1;j< rry;j++,ii++)
           for(i=rx,jj=dy >> 1;i< rrx;i++,jj++) {
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

void piramidal_decoding(int level)
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



