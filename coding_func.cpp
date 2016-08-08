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

      You should have received a copy of the GNU 
      General Public License
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

#include "def.h"
#include "globals.h"
#include "nn_search.h"
#include "prot.h"

FILE *fparam = fopen("param.txt","w");


double HV_FisherCoding(int atx, int aty, int x_size, int y_size, int *xd, int *yd, int *is, 
                                                     int *qalf, int *qbet)
{
  int x,y,i,j;
  int tip_x,tip_y,counter,isometry ;
  int isom,clas, var_class;
  int start_first, end_first, fisher_first;
  int start_second, end_second,fisher_second;
  double dist,alfa,beta, sd_r, sd_d;
  double min = 1000000000.0;
  double sum,s0;
  double mean,det;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  register double pixel;
  struct c *pointer;
  double alfa_min =0.0, beta_min=0.0, error_min = 0.0;
  tip_y =  y_size;
  tip_x =   x_size;
  int qalfa, qbeta;
     // printf("x_size = %d\t y_size = %d\n",x_size,y_size);

  // if(x_size == 1  && y_size == 1) {   /* size = 1 */
  //    // printf("x_size = %d\t y_size = %d\n",x_size,y_size);
  //    *qbet = (int)image[atx][aty];
  //    *qalf = zeroalfa;
  //    *xd = 0;
  //    *yd = 0;
  //    *is = IDENTITY;
  //    return 0.0;
  // }
  // printf("xsize= %d \t ysize = %d\n",x_size,y_size );
 
  s0 = y_size * x_size;
  for(int y=0;y < y_size;y++){
    for(int x=0;x < x_size;x++) {
       pixel = (double)image[aty+y][atx+x];
       t0 +=  pixel;
       t2 +=  pixel * pixel;
       // printf("x= %d \t y = %d\n",x,y );
       range[y][x] = pixel;
    }
  }
  // sd_r = sqrt(variance_3(x_size,y_size, range, 0,0));
  mean = t0 / s0;
      
  adaptiveNewclass_2(x_size,y_size,range,&isom,&clas); /* brings the domain in the canonical   */
  // printf("classs = %d\n",clas);  
  // printf("old Range class  = %d\n",clas);
  // flips(size,range,flip_range,isom);
  //var_class = variance_class(size,flip_range);
  // clas = getL1class(2*x_size,2*y_size,clas);
  // printf("new range class = %d\n",clas);
  // if(adaptive_fisher_class[y_size][x_size][clas] == NULL)
    // printf("adaptive_fisher_class = NULL");

  if (full_first_class) {
      start_first = 0;
      end_first   = 24;
  }
  else { 
      start_first = clas;
      end_first   = clas +1 ;
  } 

  // if (full_second_class) {
  //     start_second = 0;
  //     end_second   = 24; 
  // }
  // else { 
  //     start_second = var_class;
  //     end_second   = var_class + 1;
  // } 

  for(fisher_first = start_first; fisher_first < end_first; fisher_first++){ 
   // for(fisher_second=start_second; fisher_second < end_second; fisher_second++){ 
     pointer = HV_fisher_class[tip_y][tip_x][fisher_first]; 
     while(pointer != NULL) {
        // isometry = mapping[isom][pointer->iso];
        comparisons ++ ;
        counter ++ ;
        s1 = pointer->sum;
        s2 = pointer->sum2;
        sd_d = pointer->var;
        t1 = 0.0;
        i = pointer->ptr_x >> 1;
        j = pointer->ptr_y >> 1;
        isometry = 0;
        sum = 0.0;
        alfa = beta = 0.0;
        // printf("contract[214][216] = %u\n",contract[200][256]);
        switch(isometry) { 
         case IDENTITY   :  
                 for(x=0;x< y_size;x++)
                 for(y=0;y< x_size;y++)
                     t1 += contract[x+i][y+j]*range[x][y];
                 break;
         // case R_ROTATE90 :  
         //         for(x=0;x< y_size;x++)
         //         for(y=0;y< x_size;y++)
         //             t1 += contract[x+i][y+j]*range[y][y_size-x-1];
         //         break;
         // case L_ROTATE90 :  
         //         for(x=0;x< y_size;x++)
         //         for(y=0;y< x_size;y++)
         //             t1 += contract[x+i][y+j]*range[x_size-y-1][x];
         //         break;
         // case ROTATE180  :  
         //         for(x=0;x< y_size;x++)
         //         for(y=0;y< x_size;y++)
         //             t1 += contract[x+i][y+j]*range[y_size-x-1][x_size-y-1];
         //         break;
         // case R_VERTICAL :  
         //         for(x=0;x< y_size;x++)
         //         for(y=0;y< x_size;y++)
         //             t1 += contract[x+i][y+j]*range[x][x_size-y-1];
         //         break;
         // case R_HORIZONTAL: 
         //         for(x=0;x< y_size;x++)
         //         for(y=0;y< x_size;y++)
         //             t1 += contract[x+i][y+j]*range[y_size-x-1][y];
         //         break;
         // case F_DIAGONAL:   
         //         for(x=0;x< y_size;x++)
         //         for(y=0;y< x_size;y++)
         //             t1 += contract[x+i][y+j]*range[y][x];
         //         break;
         // case S_DIAGONAL:   
         //         for(x=0;x< y_size;x++)
         //         for(y=0;y< x_size;y++)
         //             t1 += contract[x+i][y+j]*range[x_size-y-1][y_size-x-1];
         //         break;
        }
        // printf("X size = %d \t Y size = %d\n",x_size,y_size);
        // printf("R sum = %f\n",t0);
        // printf("R sum2 = %f\n",t2);
        // printf("D sum = %f\n",s1);
        // printf("D sum2 = %f\n",s2);

        /* Compute the scalig factor */
        det = s2*s0 - s1*s1;
        if(det == 0.0)
           alfa = 0.0; 
        else
          alfa = (s0*t1 - s1*t0)  / det;
        if(alfa < 0.0 ) alfa = 0.0;
        alfa = 0.0;

        /* Quantize the scalig factor */
         qalfa =  0.5 +(alfa )/( MAX_ALFA)* (1 << N_BITALFA);
        if(qalfa < 0)  qalfa = 0;
        if(qalfa >= (1 << N_BITALFA))  qalfa = (1 << N_BITALFA) - 1;

        /* Compute the scalig factor back from the quantized value*/
        alfa = (double) qalfa / (double)(1 << N_BITALFA) * ( MAX_ALFA) ;

        /* compute offset */
        beta = (t0 - alfa * s1) / s0;    
        
        qbeta = 0.5 + beta/ ((1.0+fabs(alfa))*255)*((1<< N_BITBETA)-1);
        if (qbeta< 0) qbeta = 0;
        if (qbeta >= 1 << N_BITBETA) qbeta = (1 << N_BITBETA)-1;

        /* Compute the offset back from the quantized value*/
        beta = (double)qbeta/(double)((1 << N_BITBETA)-1)*((1.0+fabs(alfa))*255);
        if (alfa > 0.0) beta  -= alfa * 255;

        
        /* Compute the rms distance */
        sum = t2 - 2*alfa*t1 -2*beta*t0 + alfa*alfa*s2 +
                              2*alfa*beta*s1 + s0*beta*beta;
        dist = sqrt(sum / s0);

        // printf("error = %f\n",sum);
        if(dist < min ) {
            // printf("error2 = %f\n",sum);
            min = dist;
            *xd = pointer->ptr_x;
            *yd = pointer->ptr_y;
            *is = isometry;
            *qalf = qalfa;
            *qbet = qbeta;
        }
        pointer = pointer -> next;

     }
    //}
  }
   // printf("%f \t %f \t %f\n",alfa_min,beta_min,error_min);
        // printf("min=%f\n",min);

   // printf("Zero alfa count = %d\n",zero_alfa_count);
 return (min) ;
}

void HV_traverseImage(int atx, int aty, int x_size, int y_size)
{
  
  int s_size;
  int s_log;

  int max_bits = (int)log2(max_size);
  s_log = (int) log2(min_2(x_size,y_size));
  s_size = 1 << s_log;
  // printf("s_size = %d\n",s_size);

  // printf("s_log = %d\n",s_log);
  if(s_log > max_bits){
    HV_traverseImage(atx,          aty,          s_size/2, s_size/2);
    HV_traverseImage(atx+s_size/2, aty,          s_size/2, s_size/2);
    HV_traverseImage(atx,          aty+s_size/2, s_size/2, s_size/2);
    HV_traverseImage(atx+s_size/2, aty+s_size/2, s_size/2, s_size/2);
  }
  else{
    HV_compressRange(atx,aty, x_size, y_size);
  }

  if(x_size > s_size)
    HV_traverseImage(atx+s_size,aty,x_size - s_size, y_size);
  
  if(y_size > s_size)
    HV_traverseImage(atx, aty + s_size,s_size, y_size - s_size);
}

void HV_compressRange(int atx, int aty, int x_size, int y_size)
{
  // printf("isnide compress range:\n");
  printf("x_size = %d \t y_size = %d\n",x_size,y_size);
  int domx,domy,isom,qalfa,qbeta;
  static int coded = 0;
  static int oldper = 0;
  static int newper = 0;
  double compress,bpp;
  int bytes;
  int k;
  // printf("X size = %d \t Y size = %d \n",x_size,y_size);
  if(x_size == 0 || y_size == 0)
    return;

  double best_rms;
  best_rms = HV_FisherCoding(atx,aty,x_size,y_size,&domx,&domy,&isom,&qalfa,&qbeta);
  printf("Best RMS = %f\n\n",best_rms);
  if(min_2(x_size,y_size) == min_size  || best_rms < T_RMS){
    // printf("packing bits \n");
    if(max_2(x_size,y_size) !=  min_size)
      pack(1,(long)1,fp);

    // if(abs(qalfa - zeroalfa) <= zero_threshold) {
    //     qbeta = best_beta(atx,aty,x_size,0.0);
    //     qalfa = zeroalfa;
    // }
    pack(N_BITALFA, (long)qalfa, fp);
    pack(N_BITBETA, (long)qbeta, fp);
    if(isColor){
      double um = u_mean(x_size,y_size,atx,aty);
      double vm = v_mean(x_size,y_size,atx,aty);
      pack(8,(long)um,fp);
      pack(8,(long)vm,fp);
    }
    // printf("qalfa = %d\n",qalfa);
    zero_alfa_transform ++;
    if(qalfa != zeroalfa) {
       zero_alfa_transform --;
       pack(3, (long)isom, fp);
       pack(bits_per_coordinate_h,(long)(domx / SHIFT),fp);
       pack(bits_per_coordinate_w,(long)(domy / SHIFT),fp);
    }

    transforms ++;
    coded += x_size * y_size  ;
    newper = (int) (((double) coded / (image_width * image_height)) * 100);
    if(newper > oldper) {
       bytes = pack(-2,(long)0,fp);
       compress =(double) coded / (double)bytes;
       bpp      = 8.0 /compress;
       printf(" %d %% Coded, Compression %f:1, bpp %f\r", 
                                  newper,compress,bpp);
       oldper = newper;
       fflush(stdout);
    }
  }
  else{
    // printf("x size = %d \t y size = %d\n",x_size,y_size);  
    pack(1,(long)0,fp);
    int i,j,y1,y2,Rp,x1,x2,Cp;
    double hmax,vmax,v,h;
    hmax=vmax=0;
    int step = 2;
    Cp=atx;
    Rp=aty;
    for(j=aty;j<aty+y_size-1;j+=1)
    {
          y1=y2=0;
          for(i=atx;i<atx+x_size;i+=1)
          {
                y1=y1+image[j][i];
                y2=y2+image[j+1][i];
          }
          h=(min_2(j-aty+1,y_size-j-1+aty)*abs(y1-y2))/(double)(y_size);
          if(h>hmax)
          {
                hmax=h;         
                Rp=j;
          }
    }
    for(i=atx;i<atx+x_size-1;i+=1)
    {
          x1=x2=0;
          for(j=aty;j<aty+y_size;j+=1)
          {
                x1=x1+image[j][i];
                x2=x2+image[j][i+1];
          }
          v=(min_2(i+1-atx,x_size-i-1+atx)*abs(x1-x2))/(double)(x_size);
          if(v>vmax)
          {
                vmax=v;
                Cp=i;
          }
    }
    if(hmax>vmax)
    {
      pack(1,(long)1,fp);    
      pack(3,(long)Rp+1-aty,fp);    
      HV_compressRange(atx,aty,x_size,Rp-aty+1);
      HV_compressRange(atx,Rp+1,x_size,aty+y_size-Rp-1);
    }
    else
    {
      pack(1,(long)0,fp);
      pack(3,(long)Cp-atx+1,fp);
      HV_compressRange(atx,aty,Cp-atx+1,y_size);
      HV_compressRange(Cp+1,aty,atx+x_size-Cp-1,y_size);
    }
  }
}


double Saupe_FisherCoding(int atx,int aty,int size,int *xd,
                                          int *yd,int *is, int* qalf,int *qbet)
{
  int x,y,qalfa ,qbeta,i,j;
  int tip,ii,found,counter;
  int isometry, isom, clas;
  double dist,alfa,beta;
  double min = 1000000000.0;
  double sum,s0,mean;
  double det;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  static float r_vector[4096]; 
  static int nlist[MAX_NEIGHBOURS]; 
  register double pixel;


  tip = (int) rint(log((double) size)/log (2.0));

  if(tip == 0) {   /* size = 1 */
     *qbet = (int)image[atx][aty];
     *qalf = zeroalfa;
     *xd = 0;
     *yd = 0;
     *is = IDENTITY;
     return 0.0;
  }

  s0 = size * size;
  for(x=0;x < size;x++)
  for(y=0;y < size;y++) {
     pixel = (double)image[atx+x][aty+y];
     t0 +=  pixel;
     t2 +=  pixel * pixel;
     range[x][y] = pixel;
  }

  mean = t0 / s0;

  newclass(size,range,&isom,&clas);
  flips(size,range,flip_range,isom);
  ComputeSaupeVectors(flip_range,size,tip,r_vector);

  found = kdtree_search(r_vector, f_vectors[tip], feat_vect_dim[tip], 
                                        kd_tree[tip], eps ,matches,nlist); 

  for(ii=0;ii<found;ii++) {
     comparisons ++ ;
     counter ++ ;
     isometry = mapping[isom][codebook[tip][nlist[ii]].isom];
     s1 = codebook[tip][nlist[ii]].sum;
     s2 = codebook[tip][nlist[ii]].sum2;
     t1 = 0.0;
     i = codebook[tip][nlist[ii]].ptr_x >> 1;
     j = codebook[tip][nlist[ii]].ptr_y >> 1;

     switch(isometry) {
      case IDENTITY   :  
              for(x=0;x< size;x++)
              for(y=0;y< size;y++)
                  t1 += contract[x+i][y+j] * range[x][y];
              break;
      case R_ROTATE90 :  
              for(x=0;x< size;x++)
              for(y=0;y< size;y++)
                  t1 += contract[x+i][y+j] * range[y][size-x-1];
              break;
      case L_ROTATE90 :  
              for(x=0;x< size;x++)
              for(y=0;y< size;y++)
                  t1 += contract[x+i][y+j] * range[size-y-1][x];
              break;
      case ROTATE180  :  
              for(x=0;x< size;x++)
              for(y=0;y< size;y++)
                  t1 += contract[x+i][y+j] * range[size-x-1][size-y-1];
              break;
      case R_VERTICAL :  
              for(x=0;x< size;x++)
              for(y=0;y< size;y++)
                  t1 += contract[x+i][y+j] * range[x][size-y-1];
              break;
      case R_HORIZONTAL: 
              for(x=0;x< size;x++)
              for(y=0;y< size;y++)
                  t1 += contract[x+i][y+j] * range[size-x-1][y];
              break;
      case F_DIAGONAL:   
              for(x=0;x< size;x++)
              for(y=0;y< size;y++)
                  t1 += contract[x+i][y+j] * range[y][x];
              break;
      case S_DIAGONAL:   
              for(x=0;x< size;x++)
              for(y=0;y< size;y++)
                  t1 += contract[x+i][y+j] * range[size-y-1][size-x-1];
              break;
     }

     /* Compute the scalig factor */
     det = s2*s0 - s1*s1;
     if(det == 0.0)
        alfa = 0.0;
     else
       alfa = (s0*t1 - s1*t0) / det;
     if(alfa < 0.0 ) alfa = 0.0;

     /* Quantize the scalig factor */
     qalfa = 0.5 + (alfa )/( MAX_ALFA)* (1 << N_BITALFA);
     if(qalfa < 0)  qalfa = 0;
     if(qalfa >= (1 << N_BITALFA))  qalfa = (1 << N_BITALFA) - 1;

     /* Compute the scalig factor back from the quantized value*/
     alfa = (double) qalfa / (double)(1 << N_BITALFA) * ( MAX_ALFA) ;

     /* Compute the offset */
     beta= (t0 - alfa*s1) / s0;
     if (alfa > 0.0)  beta += alfa * 255;

     /* Quantize the offset */
     qbeta = 0.5 + beta/ ((1.0+fabs(alfa))*255)*((1 << N_BITBETA)-1);
     if (qbeta< 0) qbeta = 0;
     if (qbeta >= 1 << N_BITBETA) qbeta = (1 << N_BITBETA)-1;

     /* Compute the offset back from the quantized value*/
     beta = (double)qbeta/(double)((1 << N_BITBETA)-1)*((1.0+fabs(alfa))*255);
     if (alfa > 0.0) beta  -= alfa * 255;

     /* Compute the rms distance */
     sum = t2 - 2*alfa*t1 -2*beta*t0 + alfa*alfa*s2 +
                                 2*alfa*beta*s1 + s0*beta*beta;
     dist = sqrt(sum / s0);
   
     if(dist < min ) {
        min = dist;
        *xd = codebook[tip][nlist[ii]].ptr_x;
        *yd = codebook[tip][nlist[ii]].ptr_y; 
        *is = isometry;
        *qalf = qalfa;
        *qbet = qbeta;
     }
  }
 return (min) ;
}

int quan(double val)
{
  int val_int = (int)val;
  int result = 0;
  if(val - val_int < 0.5)
    result = val_int;
  else
    result = val_int + 1;

  return result;
}

double CovClass_AdaptiveSearch_FisherCoding(int atx,int aty,int size,int *xd,int *yd,int *is,
                                                       int* qalf,int *qbet)
{
  int x,y,qalfa ,qbeta,i,j;
  int tip,counter,isometry ;
  int isom,clas, stdd_class;
  int start_first, end_first, fisher_first;
  int start_second, end_second,fisher_second;
  double dist,alfa,beta, dom_var, range_var;
  double min = 1000000000.0;
  double sum,s0;
  double rmean,dmean,det;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  register double pixel;
  struct c *pointer;
  double alfa_min =0.0, beta_min=0.0, error_min = 0.0;
  tip = (int) rint((log((double) size)/log (2.0)));

  if(tip == 0) {   /* size = 1 */
     *qbet = (int)image[atx][aty];
     *qalf = zeroalfa;
     *xd = 0;
     *yd = 0;
     *is = IDENTITY;
     return 0.0;
  }

  s0 = size * size;
  for(x=0;x < size;x++)
  for(y=0;y < size;y++) {
     pixel = (double)image[atx+x][aty+y];
     t0 +=  pixel;
     t2 +=  pixel * pixel;
     range[x][y] = pixel;
  }
  range_var = variance_2(size,range,0,0);
  rmean = t0 / s0;

  newclass(size,range,&isom,&clas);
  flips(size,range,flip_range,isom);
  stdd_class = std_class(size,flip_range);
  if(stdd_class > final_max_std) stdd_class = final_max_std;

  if (full_first_class) {
      start_first = 0;
      end_first   = 3;
  }
  else { 
      start_first = clas;
      end_first   = clas + 1;
  } 

  if (full_second_class) {
      start_second = 0;
      end_second   = 24; 
  }
  else { 
      start_second = stdd_class;
      end_second   = stdd_class + 1;
  } 

  for(fisher_first = start_first; fisher_first < end_first; fisher_first ++){
    for(fisher_second=start_second; fisher_second < end_second; fisher_second++) { 
      pointer = class_std[tip][fisher_first][fisher_second]; 
      while(pointer != NULL) {
        isometry = mapping[isom][pointer->iso];
        comparisons ++ ;
        counter ++ ;
        dom_var = pointer->var;
        s1 = pointer->sum;
        s2 = pointer->sum2;
        t1 = 0.0;
        i = pointer->ptr_x >> 1;
        j = pointer->ptr_y >> 1;
        dmean = s1 / s0;
        switch(isometry) { 
         case IDENTITY   :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[x][y];
                 break;
         case R_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[y][size-x-1];
                 break;
         case L_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-y-1][x];
                 break;
         case ROTATE180  :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-x-1][size-y-1];
                 break;
         case R_VERTICAL :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[x][size-y-1];
                 break;
         case R_HORIZONTAL: 
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-x-1][y];
                 break;
         case F_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[y][x];
                 break;
         case S_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-y-1][size-x-1];
                 break;
        }

        /* Compute the scalig factor */
        det = s2*s0 - s1*s1;
        if(det == 0.0)
           alfa = 0.0;
        else
          alfa = (s0*t1 - s1*t0)  / det;
        if(alfa < 0.0 ) alfa = 0.0;
         // printf("before quantization alfa = %f\n",alfa);
        /* Quantize the scalig factor */
        // qalfa = 0.5 + (alfa )/( MAX_ALFA)* (1 << N_BITALFA);
        qalfa =  0.5 +(alfa )/( MAX_ALFA)* (1 << N_BITALFA);
        if(qalfa < 0)  qalfa = 0;
        if(qalfa >= (1 << N_BITALFA))  qalfa = (1 << N_BITALFA) - 1;

        /* Compute the scalig factor back from the quantized value*/
        alfa = (double) qalfa / (double)(1 << N_BITALFA) * ( MAX_ALFA) ;
         // printf("range_var = %f\n",range_var);            
          // printf("domain_var = %f\n",dom_var);            
        /* Condition for adaptive search */
        if((sqrt(range_var) - (alfa * sqrt(dom_var))) < (0.6 * sqrt(range_var))){    
          /* Compute the offset */
          beta= (t0 - alfa*s1) / s0;
          if (alfa > 0.0)  beta += alfa * 255;
          save_count--;
          /* Quantize the offset */
          qbeta = 0.5 + beta/ ((1.0+fabs(alfa))*255)*((1<< N_BITBETA)-1);
          if (qbeta< 0) qbeta = 0;
          if (qbeta >= 1 << N_BITBETA) qbeta = (1 << N_BITBETA)-1;
          // printf("%d\n",qbeta);
          /* Compute the offset back from the quantized value*/
          beta = (double)qbeta/(double)((1 << N_BITBETA)-1)*((1.0+fabs(alfa))*255);
          if (alfa > 0.0) beta  -= alfa * 255;

          // dist = sqrt(range_var + alfa * alfa * dom_var - 
                       // 2 * alfa * (rmean * dmean - (t1/s0)));
          /* Compute the rms distance */
          sum = t2 - 2*alfa*t1 -2*beta*t0 + alfa*alfa*s2 +
                                2*alfa*beta*s1 + s0*beta*beta;
          // sum = range_var + alfa * alfa * dom_var - (2 * alfa * ((t0/s0 * s1/s0) - t1/(s0 * s0))); 

          dist = sqrt(sum / s0);
           
          if(dist < min ) {
              // if(alfa == 0.0) 
                // printf("alfa = %f\n",alfa);
              
              min = dist;
              // error_min = dist;
              // alfa_min = alfa;
              // beta_min = beta;
              *xd = pointer->ptr_x;
              *yd = pointer->ptr_y;
              *is = isometry;
              *qalf = qalfa;
              *qbet = qbeta;
          }
        }
        pointer = pointer -> next;
      }
    }
  }
 // printf("%f \t %f \t %f\n",alfa_min,beta_min,error_min);
 // printf("Zero alfa count = %d\n",zero_alfa_count);
 return (min) ;
}

double AdaptiveSearch_FisherCoding(int atx,int aty,int size,int *xd,int *yd,int *is,
                                                       int* qalf,int *qbet)
{
  int x,y,qalfa ,qbeta,i,j;
  int tip,counter,isometry ;
  int isom,clas, var_class;
  int start_first, end_first, fisher_first;
  int start_second, end_second,fisher_second;
  double dist,alfa,beta, dom_var, range_var;
  double min = 1000000000.0;
  double sum,s0;
  double mean,det;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  register double pixel;
  struct c *pointer;
  double alfa_min =0.0, beta_min=0.0, error_min = 0.0;
  tip = (int) rint((log((double) size)/log (2.0)));

  if(tip == 0) {   /* size = 1 */
     *qbet = (int)image[atx][aty];
     *qalf = zeroalfa;
     *xd = 0;
     *yd = 0;
     *is = IDENTITY;
     return 0.0;
  }

  s0 = size * size;
  for(x=0;x < size;x++)
  for(y=0;y < size;y++) {
     pixel = (double)image[atx+x][aty+y];
     t0 +=  pixel;
     t2 +=  pixel * pixel;
     range[x][y] = pixel;
  }
  range_var = variance(size,size,atx,aty);
  mean = t0 / s0;

  newclass(size,range,&isom,&clas);
  flips(size,range,flip_range,isom);
  var_class = variance_class(size,flip_range);

  if (full_first_class) {
      start_first = 0;
      end_first   = 3;
  }
  else { 
      start_first = clas;
      end_first   = clas + 1;
  } 

  if (full_second_class) {
      start_second = 0;
      end_second   = 24; 
  }
  else { 
      start_second = var_class;
      end_second   = var_class + 1;
  } 

  for(fisher_first = start_first; fisher_first < end_first; fisher_first ++){
    for(fisher_second=start_second; fisher_second < end_second; fisher_second++) { 
      pointer = class_fisher[tip][fisher_first][fisher_second]; 
      while(pointer != NULL) {
        isometry = mapping[isom][pointer->iso];
        comparisons ++ ;
        counter ++ ;
        dom_var = pointer->var;
        s1 = pointer->sum;
        s2 = pointer->sum2;
        t1 = 0.0;
        i = pointer->ptr_x >> 1;
        j = pointer->ptr_y >> 1;
   
        switch(isometry) { 
         case IDENTITY   :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[x][y];
                 break;
         case R_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[y][size-x-1];
                 break;
         case L_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-y-1][x];
                 break;
         case ROTATE180  :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-x-1][size-y-1];
                 break;
         case R_VERTICAL :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[x][size-y-1];
                 break;
         case R_HORIZONTAL: 
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-x-1][y];
                 break;
         case F_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[y][x];
                 break;
         case S_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-y-1][size-x-1];
                 break;
        }

        /* Compute the scalig factor */
        det = s2*s0 - s1*s1;
        if(det == 0.0)
           alfa = 0.0;
        else
          alfa = (s0*t1 - s1*t0)  / det;
        if(alfa < 0.0 ) alfa = 0.0;
         // printf("before quantization alfa = %f\n",alfa);
        /* Quantize the scalig factor */
        // qalfa = 0.5 + (alfa )/( MAX_ALFA)* (1 << N_BITALFA);
        qalfa =  0.5 +(alfa )/( MAX_ALFA)* (1 << N_BITALFA);
        if(qalfa < 0)  qalfa = 0;
        if(qalfa >= (1 << N_BITALFA))  qalfa = (1 << N_BITALFA) - 1;

        /* Compute the scalig factor back from the quantized value*/
        alfa = (double) qalfa / (double)(1 << N_BITALFA) * ( MAX_ALFA) ;
        
        /* Condition for adaptive search */
        if(sqrt(range_var) - alfa * sqrt(dom_var) < 0.6 * sqrt(range_var)){    
          
          /* Compute the offset */
          beta= (t0 - alfa*s1) / s0;
          if (alfa > 0.0)  beta += alfa * 255;

          /* Quantize the offset */
          qbeta = 0.5 + beta/ ((1.0+fabs(alfa))*255)*((1<< N_BITBETA)-1);
          if (qbeta< 0) qbeta = 0;
          if (qbeta >= 1 << N_BITBETA) qbeta = (1 << N_BITBETA)-1;
          // printf("%d\n",qbeta);
          /* Compute the offset back from the quantized value*/
          beta = (double)qbeta/(double)((1 << N_BITBETA)-1)*((1.0+fabs(alfa))*255);
          if (alfa > 0.0) beta  -= alfa * 255;

          /* Compute the rms distance */
          sum = t2 - 2*alfa*t1 -2*beta*t0 + alfa*alfa*s2 +
                                2*alfa*beta*s1 + s0*beta*beta;
          // sum = range_var + alfa * alfa * dom_var - (2 * alfa * ((t0/s0 * s1/s0) - t1/(s0 * s0))); 

          dist = sqrt(sum / s0);
           
          if(dist < min ) {
              min = dist;
              *xd = pointer->ptr_x;
              *yd = pointer->ptr_y;
              *is = isometry;
              *qalf = qalfa;
              *qbet = qbeta;
          }
        }
        pointer = pointer -> next;
      }
    }
  }
 // printf("%f \t %f \t %f\n",alfa_min,beta_min,error_min);
 // printf("Zero alfa count = %d\n",zero_alfa_count);
 return (min) ;
}

double testing_FisherCoding(int atx,int aty,int size,int *xd,int *yd,int *is,
                                                       int* qalf,int *qrmea)
{
  int x,y,qalfa ,qbeta,i,j,qrmean;
  int tip,counter,isometry ;
  int isom,clas, var_class;
  int start_first, end_first, fisher_first;
  int start_second, end_second,fisher_second;
  double dist,alfa,beta, rmean;
  double min = 1000000000.0;
  double sum,s0;
  double mean,det;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  register double pixel;
  struct c *pointer;
  double alfa_min =0.0, beta_min=0.0, error_min = 0.0;
  tip = (int) rint((log((double) size)/log (2.0)));

  if(tip == 0) {   /* size = 1 */
     *qrmea = (int)image[atx][aty];
     *qalf = zeroalfa;
     *xd = 0;
     *yd = 0;
     *is = IDENTITY;
     return 0.0;
  }

  s0 = size * size;
  for(x=0;x < size;x++)
  for(y=0;y < size;y++) {
     pixel = (double)image[atx+x][aty+y];
     t0 +=  pixel;
     t2 +=  pixel * pixel;
     range[x][y] = pixel;
  }

  mean = t0 / s0;

  newclass(size,range,&isom,&clas);
  flips(size,range,flip_range,isom);
  var_class = variance_class(size,flip_range);

  if (full_first_class) {
      start_first = 0;
      end_first   = 3;
  }
  else { 
      start_first = clas;
      end_first   = clas + 1;
  } 

  if (full_second_class) {
      start_second = 0;
      end_second   = 24; 
  }
  else { 
      start_second = var_class;
      end_second   = var_class + 1;
  } 

  for(fisher_first = start_first; fisher_first < end_first; fisher_first ++) 
  for(fisher_second=start_second; fisher_second < end_second; fisher_second++) { 
     pointer = class_fisher[tip][fisher_first][fisher_second]; 
     while(pointer != NULL) {
        isometry = mapping[isom][pointer->iso];
        comparisons ++ ;
        counter ++ ;
        s1 = pointer->sum;
        s2 = pointer->sum2;
       
        t1 = 0.0;
        i = pointer->ptr_x >> 1;
        j = pointer->ptr_y >> 1;
   
        switch(isometry) { 
         case IDENTITY   :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[x][y];
                 break;
         case R_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[y][size-x-1];
                 break;
         case L_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-y-1][x];
                 break;
         case ROTATE180  :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-x-1][size-y-1];
                 break;
         case R_VERTICAL :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[x][size-y-1];
                 break;
         case R_HORIZONTAL: 
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-x-1][y];
                 break;
         case F_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[y][x];
                 break;
         case S_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-y-1][size-x-1];
                 break;
        }

        /* Compute the scalig factor */
        det = s2*s0 - s1*s1;
        if(det == 0.0)
           alfa = 0.0;
        else
          alfa = (s0*t1 - s1*t0)  / det;
        if(alfa < 0.0 ) alfa = 0.0;
         // printf("before quantization alfa = %f\n",alfa);
        /* Quantize the scalig factor */
        // qalfa = 0.5 + (alfa )/( MAX_ALFA)* (1 << N_BITALFA);
        qalfa =  0.5 +(alfa )/( MAX_ALFA)* (1 << N_BITALFA);
        if(qalfa < 0)  qalfa = 0;
        if(qalfa >= (1 << N_BITALFA))  qalfa = (1 << N_BITALFA) - 1;

        /* Compute the scalig factor back from the quantized value*/
        alfa = (double) qalfa / (double)(1 << N_BITALFA) * ( MAX_ALFA) ;
//      
        if(alfa == 0)
          zero_alfa_count++;
        // printf("after deqauantization alfa = %f\n\n",alfa);

        /* Compute the offset */
        beta= (t0 - alfa*s1) / s0;
        if (alfa > 0.0)  beta += alfa * 255;

        /* Quantize the offset */
        qbeta = 0.5 + beta/ ((1.0+fabs(alfa))*255)*((1<< N_BITBETA)-1);
        if (qbeta< 0) qbeta = 0;
        if (qbeta >= 1 << N_BITBETA) qbeta = (1 << N_BITBETA)-1;
        // printf("%d\n",qbeta);
        /* Compute the offset back from the quantized value*/
        beta = (double)qbeta/(double)((1 << N_BITBETA)-1)*((1.0+fabs(alfa))*255);
        if (alfa > 0.0) beta  -= alfa * 255;
        
        rmean = beta + alfa*s1/s0;
        qrmean = 0.5 + rmean/ ((1.0+fabs(alfa))*255)*((1<< N_BITRMEAN)-1);
        if (qrmean >= 1 << N_BITRMEAN) qrmean = (1 << N_BITRMEAN)-1;
        /* Compute the rms distance */
       // sum = t2 - 2*alfa*t1 -2*beta*t0 + alfa*alfa*s2 +
                           //   2*alfa*beta*s1 + s0*beta*beta;  //original error
        // sum = range_var + alfa * alfa * dom_var - (2 * alfa * ((t0/s0 * s1/s0) - t1/(s0 * s0))); 
        sum = t2 + (s0 * t0 * t0) - 2*(t0 * t0) - (alfa * s2/s0) - (alfa * s1/s0 * s1/s0)*s0 + 2*(alfa * s1/s0 * s1/s0)
               - 2 *(t1/s0 * alfa) + 4*(alfa * t0/s0 * s1/s0) - 2*(alfa * t0/s0 * s1);  
        dist = sqrt(sum);
         
        if(dist < min ) {
            min = dist;
            error_min = dist;
            alfa_min = alfa;
            beta_min = beta;
            *xd = pointer->ptr_x;
            *yd = pointer->ptr_y;
            *is = isometry;
            *qalf = qalfa;
            *qrmea = qrmean;
            //*qrmea = qbeta;
        }
        pointer = pointer -> next;
        //printf("%f\n",error_min);

     }
  }

// printf("%f \t %f \t %f\n",alfa_min,beta_min,error_min);
  // printf("Zero alfa count = %d\n",zero_alfa_count);
 return (min) ;
}

double LumInv_FisherCoding(int atx,int aty,int size,int *xd,int *yd,int *is,
                                                       int* qalf,int* qrmea)
{
  int x,y,qalfa ,qrmean,i,j,qbeta;
  int tip,counter,isometry ;
  int isom,clas, var_class;
  int start_first, end_first, fisher_first;
  int start_second, end_second,fisher_second;
  double dist,alfa, rmean_min, rmean,beta;
  double min = 1000000000.0;
  double sum,s0;
  double mean,det;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  register double pixel;
  struct c *pointer;
  double alfa_min =0.0, beta_min=0.0, error_min = 0.0;
  tip = (int) rint((log((double) size)/log (2.0)));

  if(tip == 0) {   /* size = 1 */
     *qrmea = (int)image[atx][aty];
     *qalf = zeroalfa;
     *xd = 0;
     *yd = 0;
     *is = IDENTITY;
     return 0.0;
  }

  s0 = size * size;
  for(x=0;x < size;x++)
  for(y=0;y < size;y++) {
     pixel = (double)image[atx+x][aty+y];
     t0 +=  pixel;
     t2 +=  pixel * pixel;
     range[x][y] = pixel;
  }

  mean = t0 / s0;

  newclass(size,range,&isom,&clas);
  flips(size,range,flip_range,isom);
  var_class = variance_class(size,flip_range);

  if (full_first_class) {
      start_first = 0;
      end_first   = 3;
  }
  else { 
      start_first = clas;
      end_first   = clas + 1;
  } 

  if (full_second_class) {
      start_second = 0;
      end_second   = 24; 
  }
  else { 
      start_second = var_class;
      end_second   = var_class + 1;
  } 

  for(fisher_first = start_first; fisher_first < end_first; fisher_first ++) 
  for(fisher_second=start_second; fisher_second < end_second; fisher_second++) { 
     pointer = class_fisher[tip][fisher_first][fisher_second]; 
     while(pointer != NULL) {
        isometry = mapping[isom][pointer->iso];
        comparisons ++ ;
        counter ++ ;
        s1 = pointer->sum;
        s2 = pointer->sum2;
       
        t1 = 0.0;
        i = pointer->ptr_x >> 1;
        j = pointer->ptr_y >> 1;
   
        switch(isometry) { 
         case IDENTITY   :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[x][y];
                 break;
         case R_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[y][size-x-1];
                 break;
         case L_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-y-1][x];
                 break;
         case ROTATE180  :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-x-1][size-y-1];
                 break;
         case R_VERTICAL :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[x][size-y-1];
                 break;
         case R_HORIZONTAL: 
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-x-1][y];
                 break;
         case F_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[y][x];
                 break;
         case S_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-y-1][size-x-1];
                 break;
        }

        /* Compute the scalig factor */
        det = s2*s0 - s1*s1;
        if(det == 0.0)
           alfa = 0.0;
        else
          alfa = (s0*t1 - s1*t0) / det;
        if(alfa < 0.0 ) alfa = 0.0;
       
        qalfa =  0.5 +(alfa )/( MAX_ALFA)* (1 << N_BITALFA);
        if(qalfa < 0)  qalfa = 0;
        if(qalfa >= (1 << N_BITALFA))  qalfa = (1 << N_BITALFA) - 1;

        /* Compute the scalig factor back from the quantized value*/
        alfa = (double) qalfa / (double)(1 << N_BITALFA) * ( MAX_ALFA) ;
        // if(alfa == 0)
          // zero_alfa_count++;
        // printf("deq alfa = %f\n",alfa);
        /*quantize the range mean */
       
        
      
        beta= (t0 - alfa*s1) / s0;
        if (alfa > 0.0)  beta += alfa * 255;

        /* Quantize the offset */
        qbeta = 0.5 + beta/ ((1.0+fabs(alfa))*255)*((1<< N_BITBETA)-1);
        if (qbeta< 0) qbeta = 0;
        if (qbeta >= 1 << N_BITBETA) qbeta = (1 << N_BITBETA)-1;

        /* Compute the offset back from the quantized value*/
        beta = (double)qbeta/(double)((1 << N_BITBETA)-1)*((1.0+fabs(alfa))*255);
        if (alfa > 0.0) beta  -= alfa * 255;
        /* Compute the rms distance */
        sum = t2 - 2*alfa*t1 -2*beta*t0 + alfa*alfa*s2 +
                              2*alfa*beta*s1 + s0*beta*beta;
        dist = sqrt(sum / s0);
         printf("%d\n",qbeta);
        if(dist < min ) {
            min = dist;
            error_min = dist;
            rmean = beta + (alfa * s1)/s0;
            qrmean = quan(rmean);
            alfa_min = alfa;
            rmean_min = rmean;
            *xd = pointer->ptr_x;
            *yd = pointer->ptr_y;
            *is = isometry;
            *qalf = qalfa;
            *qrmea = qbeta;
        }
        pointer = pointer -> next;
     }
  }

// printf("%f \t %f \t %f\n",alfa_min,rmean_min,Ferror_min);
  // printf("zero alfa count = %d\n",zero_alfa_count);
 return (min) ;
}

double Nonlinear_FisherCoding(int atx,int aty,int size,int *xd,int *yd,int *is,
                                                       int* qalf1,int* qalf2,int *qbet)
{
  int x,y,qalfa1, qalfa2 ,qbeta,i,j;
  int tip,counter,isometry ;
  int isom,clas, var_class;
  int start_first, end_first, fisher_first;
  int start_second, end_second,fisher_second;
  double dist,alfa1,alfa2,beta;
  double min = 1000000000.0;
  double sum,s0;
  double mean,det1, det2;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double t3 = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  double s3 = 0.0;
  double s4 = 0.0;
  double test_s1 = 0.0;
double alfa1_min =0.0, beta_min=0.0, error_min = 0.0,alfa2_min = 0.0;

  register double pixel;
  struct c *pointer;

  tip = (int) rint((log((double) size)/log (2.0)));

  if(tip == 0) {   /* size = 1 */
     *qbet = (int)image[atx][aty];
     *qalf1 = zeroalfa;
     *qalf2 = zeroalfa;
     *xd = 0;
     *yd = 0;
     *is = IDENTITY;
     return 0.0;
  }
  
  s0 = size * size;
  for(x=0;x < size;x++)
  for(y=0;y < size;y++) {
     pixel = (double)image[atx+x][aty+y];
     t0 +=  pixel;
     t2 +=  pixel * pixel;
     range[x][y] = pixel;
  }

  mean = t0 / s0;
  //printf("size= %d\n",size);
  newclass(size,range,&isom,&clas);
  flips(size,range,flip_range,isom);
  var_class = variance_class(size,flip_range);

  if (full_first_class) {
      start_first = 0;
      end_first   = 3;
  }
  else { 
      start_first = clas;
      end_first   = clas + 1;
  } 

  if (full_second_class) {
      start_second = 0;
      end_second   = 24; 
  }
  else { 
      start_second = var_class;
      end_second   = var_class + 1;
  } 

  for(fisher_first = start_first; fisher_first < end_first; fisher_first ++) 
  for(fisher_second=start_second; fisher_second < end_second; fisher_second++) { 
     pointer = class_fisher[tip][fisher_first][fisher_second]; 
     while(pointer != NULL) {
        isometry = mapping[isom][pointer->iso];
        comparisons ++ ;
        counter ++ ;
        s1 = pointer->sum ;
        s2 = pointer->sum2;
        s3 = pointer->sum3;
        s4 = pointer->sum4;
      
        t1 = 0.0;
        t3 = 0.0;
        test_s1 = 0.0;
        i = pointer->ptr_x >> 1;
        j = pointer->ptr_y >> 1;
   
        switch(isometry) { 
         case IDENTITY   :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++){
                     t1 += contract[x+i][y+j]*range[x][y];
                     test_s1 += contract[x+i][y+j];
                     t3 += contract[x+i][y+j]*contract[x+i][y+j]*range[x][y];
                  }
                 break;
         case R_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++){
                     t1 += contract[x+i][y+j]*range[y][size-x-1];
                     test_s1 += contract[x+i][y+j];
                     t3 += contract[x+i][y+j]*contract[x+i][y+j]*range[y][size-x-1];
                  }

                 break;
         case L_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++){
                     t1 += contract[x+i][y+j]*range[size-y-1][x];
                     test_s1 += contract[x+i][y+j];
                     t3 += contract[x+i][y+j]*range[size-y-1][x]*contract[x+i][y+j];
                  }
                 break;
         case ROTATE180  :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++){
                    t1 += contract[x+i][y+j]*range[size-x-1][size-y-1];
                    test_s1 += contract[x+i][y+j];
                    t3 += contract[x+i][y+j]*range[size-x-1][size-y-1]*contract[x+i][y+j];
                  }  
                 break;
         case R_VERTICAL :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++){
                     t1 += contract[x+i][y+j]*range[x][size-y-1];
                     test_s1 += contract[x+i][y+j];
                     t3 += contract[x+i][y+j]*range[x][size-y-1]*contract[x+i][y+j];
                  }
                 break;
         case R_HORIZONTAL: 
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++){
                     t1 += contract[x+i][y+j]*range[size-x-1][y];
                     test_s1 += contract[x+i][y+j];
                     t3 += contract[x+i][y+j]*range[size-x-1][y]*contract[x+i][y+j];
                  }
                 break;
         case F_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++){
                     t1 += contract[x+i][y+j]*range[y][x];
                     test_s1 += contract[x+i][y+j];
                     t3 += contract[x+i][y+j]*range[y][x]*contract[x+i][y+j];
                  }
                 break;
         case S_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++){
                     t1 += contract[x+i][y+j]*range[size-y-1][size-x-1];
                     test_s1 += contract[x+i][y+j];
                     t3 += contract[x+i][y+j]*range[size-y-1][size-x-1]*contract[x+i][y+j];
                  }
                 break;
        }

        /* Compute the scalig factor */
      
        double r = t0/s0;
        double d = s1/s0;
        double d2 = s2/s0;
        double d3 = s3/s0;
        double d4 = s4/s0;
        double rd = t1/s0;
        double rd2 = t3/s0;
        int flag1 = 0, flag2 = 0;
        
        det1 = (((d2*d) - d3) * ((d2*d) - d3)) - (((d*d) - d2) * ((d2 * d2) - (d*d*d*d))) ;
        
        if(det1 == 0.0)
          alfa1 = 0.0;
        else
          alfa1 = 1000*((rd - r*d) * (d2 * d2 - d4) - (rd2 - d2*r)*(d2*d - d3))/det1;
        if(alfa1 < 0) alfa1   = 0;
           
         // printf("alfa1 = %f\n",alfa1);

        qalfa1 = 0.5 + (alfa1 )/( MAX_ALFA1)* (1 << N_BITALFA1);
        if(qalfa1 < 0)  qalfa1 = 0;
        if(qalfa1 >= (1 << N_BITALFA1))  qalfa1 = (1 << N_BITALFA1) - 1;
        

        alfa1 = (double) qalfa1 / (double)(1 << N_BITALFA1) * (MAX_ALFA1) ;
         // printf("deq alfa1 = %f\n",alfa1);
        
          
        det2 = d4 - (d2 * d2);
       
        if(det2 == 0.0)
          alfa2 = 0.0;
        else
          alfa2 = ((rd2 - r*d2) + alfa1*(d*d2 - d3))/det2 ;
        if(alfa2 < 0.0) alfa2  = 0;
         // printf("alfa2 = %f\n",alfa2);
       
       qalfa2 = 0.005 + (alfa2 )/( MAX_ALFA2)* (1 << N_BITALFA2);
        if(qalfa2 < 0)  qalfa2 = 0;
        if(qalfa2 >= (1 << N_BITALFA2))  qalfa2 = (1 << N_BITALFA2) - 1;
        
        alfa2 = (double) qalfa2 / (double)(1 << N_BITALFA2) * ( MAX_ALFA2) ;
         // printf("deq alfa2 = %f\n\n",alfa2);
        
    
    //    int flag_beta = 0;
        beta = r - (alfa1 * d) - (alfa2 * d2);
        if(beta <0.0 ) beta = 0.0;
        if (alfa1 > 0.0)  beta += alfa1 * 255 ;
       
      
        qbeta = 0.5 + beta/ ((1.0+fabs(alfa1))*255)*((1<< N_BITBETA2)-1);
        if (qbeta< 0) qbeta = 0;
        if (qbeta >= 1 << N_BITBETA2) qbeta = (1 << N_BITBETA2)-1;
    

        beta = (double)qbeta/(double)((1 << N_BITBETA2)-1)*((1.0+fabs(alfa1))*255);
        if (alfa1 > 0.0) beta  -= alfa1 * 255;


       sum = t2 + (s0 * beta) - (2 * beta * t0) + (alfa1 * alfa1 * s2) + (alfa2 * alfa2 * s4)
               + (alfa1 * alfa2 * s3) - (2 * alfa1 * t1) - (2 * alfa2 * t3)
               + (2 * alfa1 * beta * s1) + (2 * alfa2 * beta * s2);   


       dist = sqrt(sum/s0);
        
       //
        if(dist < min ) {
            min = dist;
            error_min = min;
            alfa1_min = alfa1;
            alfa2_min = alfa2;
            beta_min = beta;
            *xd = pointer->ptr_x;
            *yd = pointer->ptr_y;
            *is = isometry;
            *qalf1 = qalfa1;
            *qalf2 = qalfa2;
            *qbet = qbeta;
        }
        pointer = pointer -> next;
     }
  }
 // printf("%f \t %f \t %f \t %f\n",alfa1_min,alfa2_min,beta_min,error_min);
 // printf("error: %f\n",error_min);
 return (min) ;
}

double EntropyCoding(int atx,int aty,int size,int *xd,int *yd,int *is,
                                                       int* qalf,int *qbet)
{
  int x,y,qalfa ,qbeta,i,j;
  int tip,counter,isometry ;
  int isom,clas, entr_class;
  int start_first, end_first, fisher_first;
  int start_second, end_second,fisher_second;
  double dist,alfa,beta;
  double min = 1000000000.0;
  double sum,s0;
  double mean,det;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  register double pixel;
  struct c *pointer;
  double alfa_min =0.0, beta_min=0.0, error_min = 0.0;
  tip = (int) rint((log((double) size)/log (2.0)));

  if(tip == 0) {   /* size = 1 */
     *qbet = (int)image[atx][aty];
     *qalf = zeroalfa;
     *xd = 0;
     *yd = 0;
     *is = IDENTITY;
     return 0.0;
  }

  s0 = size * size;
  for(x=0;x < size;x++)
  for(y=0;y < size;y++) {
     pixel = (double)image[atx+x][aty+y];
     t0 +=  pixel;
     t2 +=  pixel * pixel;
     range[x][y] = pixel;
  }
  mean = t0 / s0;

  newclass(size,range,&isom,&clas);
  flips(size,range,flip_range,isom);
  entr_class = ent_class(size,flip_range);
  if (full_first_class) {
      start_first = 0;
      end_first   = 3;
  }
  else { 
      start_first = clas;
      end_first   = clas + 1;
  } 

  if (full_second_class) {
      start_second = 0;
      end_second   = 24; 
  }
  else { 
      start_second = entr_class;
      end_second   = entr_class + 1;
  } 

  for(fisher_first = start_first; fisher_first < end_first; fisher_first ++) 
  for(fisher_second=start_second; fisher_second < end_second; fisher_second++) { 
     pointer = class_entropy[tip][fisher_first][fisher_second];
     // printf("tip = %d, fisher_first=%d, fisher_second=%d\n",tip,fisher_first,fisher_second); 
     while(pointer != NULL) {
      // printf("inside while\n");
        isometry = mapping[isom][pointer->iso];
        comparisons ++ ;
        counter ++ ;
        s1 = pointer->sum;
        s2 = pointer->sum2;
        t1 = 0.0;
        i = pointer->ptr_x >> 1;
        j = pointer->ptr_y >> 1;
   
        switch(isometry) { 
         case IDENTITY   :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += range[x][y];
                 break;
         case R_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[y][size-x-1];
                 break;
         case L_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-y-1][x];
                 break;
         case ROTATE180  :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-x-1][size-y-1];
                 break;
         case R_VERTICAL :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[x][size-y-1];
                 break;
         case R_HORIZONTAL: 
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-x-1][y];
                 break;
         case F_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[y][x];
                 break;
         case S_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-y-1][size-x-1];
                 break;
        }

        /* Compute the scalig factor */
        det = s2*s0 - s1*s1;
        if(det == 0.0)
           alfa = 0.0;
        else
          alfa = (s0*t1 - s1*t0)  / det;
        if(alfa < 0.0 ) alfa = 0.0;

        /* Quantize the scalig factor */
        qalfa =  0.5 +(alfa )/( MAX_ALFA)* (1 << N_BITALFA);
        if(qalfa < 0)  qalfa = 0;
        if(qalfa >= (1 << N_BITALFA))  qalfa = (1 << N_BITALFA) - 1;

        /* Compute the scalig factor back from the quantized value*/
        alfa = (double) qalfa / (double)(1 << N_BITALFA) * ( MAX_ALFA) ;
        // printf("after deqauantization alfa = %f\n\n",alfa);

        /* Compute the offset */
        beta= (t0 - alfa*s1) / s0;
        if (alfa > 0.0)  beta += alfa * 255;

        /* Quantize the offset */
        qbeta = 0.5 + beta/ ((1.0+fabs(alfa))*255)*((1<< N_BITBETA)-1);
        if (qbeta< 0) qbeta = 0;
        if (qbeta >= 1 << N_BITBETA) qbeta = (1 << N_BITBETA)-1;
        // printf("%d\n",qbeta);
        /* Compute the offset back from the quantized value*/
        beta = (double)qbeta/(double)((1 << N_BITBETA)-1)*((1.0+fabs(alfa))*255);
        if (alfa > 0.0) beta  -= alfa * 255;

        /* Compute the rms distance */
        sum = t2 - 2*alfa*t1 -2*beta*t0 + alfa*alfa*s2 +
                              2*alfa*beta*s1 + s0*beta*beta;
        dist = sqrt(sum / s0);
         
        if(dist < min ) {
            min = dist;
            error_min = dist;
            alfa_min = alfa;
            beta_min = beta;
            *xd = pointer->ptr_x;
            *yd = pointer->ptr_y;
            *is = isometry;
            *qalf = qalfa;
            *qbet = qbeta;
        }
        pointer = pointer -> next;
        // printf("%f\n",error_min);

     }
  }

// printf("%f \t %f \t %f\n",alfa_min,beta_min,error_min);
  // printf("Zero alfa count = %d\n",zero_alfa_count);
 return (min) ;
}

double COVCoding(int atx,int aty,int size,int *xd,int *yd,int *is,
                                                       int* qalf,int *qbet)
{
  int x,y,qalfa ,qbeta,i,j;
  int tip,counter,isometry ;
  int isom,clas, stdd_class;
  int start_first, end_first, fisher_first;
  int start_second, end_second,fisher_second;
  double dist,alfa,beta, sd_r, sd_d;
  double min = 1000000000.0;
  double sum,s0;
  double mean,det;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  register double pixel;
  struct c *pointer;
  double alfa_min =0.0, beta_min=0.0, error_min = 0.0;
  tip = (int) rint((log((double) size)/log (2.0)));

  if(tip == 0) {   /* size = 1 */
     *qbet = (int)image[atx][aty];
     *qalf = zeroalfa;
     *xd = 0;
     *yd = 0;
     *is = IDENTITY;
     return 0.0;
  }

  s0 = size * size;
  for(x=0;x < size;x++)
  for(y=0;y < size;y++) {
     pixel = (double)image[atx+x][aty+y];
     t0 +=  pixel;
     t2 +=  pixel * pixel;
     range[x][y] = pixel;
  }
  mean = t0 / s0;

  newclass(size,range,&isom,&clas);
  flips(size,range,flip_range,isom);
  stdd_class = std_class(size,flip_range);
  if(stdd_class > final_max_std) stdd_class = final_max_std;
  
  if (full_first_class) {
      start_first = 0;
      end_first   = 3;
  }
  else { 
      start_first = clas;
      end_first   = clas + 1;
  } 

  if (full_second_class) {
      start_second = 0;
      end_second   = 24; 
  }
  else { 
      start_second = stdd_class;
      end_second   = stdd_class + 1;
  } 

  for(fisher_first = start_first; fisher_first < end_first; fisher_first ++) 
  for(fisher_second=start_second; fisher_second < end_second; fisher_second++) { 
     pointer = class_std[tip][fisher_first][fisher_second];
     // printf("tip = %d, fisher_first=%d, fisher_second=%d\n",tip,fisher_first,fisher_second); 
     while(pointer != NULL) {
      // printf("inside while\n");
      // printf("isom : %d \t pointer->iso : %d\n",isom,pointer->iso);
      if(pointer->iso != NULL)
        isometry = mapping[isom][pointer->iso];
        comparisons ++ ;
        counter ++ ;
        s1 = pointer->sum;
        s2 = pointer->sum2;
        sd_d = pointer->var;
        t1 = 0.0;
        i = pointer->ptr_x >> 1;
        j = pointer->ptr_y >> 1;
   
        switch(isometry) { 
         case IDENTITY   :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[x][y];
                 break;
         case R_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[y][size-x-1];
                 break;
         case L_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-y-1][x];
                 break;
         case ROTATE180  :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-x-1][size-y-1];
                 break;
         case R_VERTICAL :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[x][size-y-1];
                 break;
         case R_HORIZONTAL: 
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-x-1][y];
                 break;
         case F_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[y][x];
                 break;
         case S_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-y-1][size-x-1];
                 break;
        }

        /* Compute the scalig factor */
        det = s2*s0 - s1*s1;  
        if(det == 0.0)
           alfa = 0.0;
        else
          alfa = (s0*t1 - s1*t0)  / det;
        if(alfa < 0.0 ) alfa = 0.0;

        /* Quantize the scalig factor */
        qalfa =  0.5 +(alfa )/( MAX_ALFA)* (1 << N_BITALFA);
        if(qalfa < 0)  qalfa = 0;
        if(qalfa >= (1 << N_BITALFA))  qalfa = (1 << N_BITALFA) - 1;

        /* Compute the scalig factor back from the quantized value*/
        alfa = (double) qalfa / (double)(1 << N_BITALFA) * ( MAX_ALFA) ;
        // printf("after deqauantization alfa = %f\n\n",alfa);

        /* Compute the offset */
        beta= (t0 - alfa*s1) / s0;
        if (alfa > 0.0)  beta += alfa * 255;

        /* Quantize the offset */
        qbeta = 0.5 + beta/ ((1.0+fabs(alfa))*255)*((1<< N_BITBETA)-1);
        if (qbeta< 0) qbeta = 0;
        if (qbeta >= 1 << N_BITBETA) qbeta = (1 << N_BITBETA)-1;
        // printf("%d\n",qbeta);
        /* Compute the offset back from the quantized value*/
        beta = (double)qbeta/(double)((1 << N_BITBETA)-1)*((1.0+fabs(alfa))*255);
        if (alfa > 0.0) beta  -= alfa * 255;

        /* Compute the rms distance */
        sum = t2 - 2*alfa*t1 -2*beta*t0 + alfa*alfa*s2 +
                              2*alfa*beta*s1 + s0*beta*beta;
        dist = sqrt(sum / s0);
         
        if(dist < min ) {
            min = dist;
            error_min = dist;
            alfa_min = alfa;
            beta_min = beta;
            *xd = pointer->ptr_x;
            *yd = pointer->ptr_y;
            *is = isometry;
            *qalf = qalfa;
            *qbet = qbeta;
        }
        pointer = pointer -> next;
        // printf("%f\nN",error_min);

     }
  }

// printf("%f \t %f \t %f\n",alfa_min,beta_min,error_min);
  // printf("Zero alfa count = %d\n",zero_alfa_count);
 return (min) ;
}

double STDCoding(int atx,int aty,int size,int *xd,int *yd,int *is,
                                                       int* qalf,int *qbet)
{
  int x,y,qalfa ,qbeta,i,j;
  int tip,counter,isometry ;
  int isom,clas, stdd_class;
  int start_first, end_first, fisher_first;
  int start_second, end_second,fisher_second;
  double dist,alfa,beta, sd_r, sd_d;
  double min = 1000000000.0;
  double sum,s0;
  double mean,det;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  register double pixel;
  struct c *pointer;
  double alfa_min =0.0, beta_min=0.0, error_min = 0.0;
  tip = (int) rint((log((double) size)/log (2.0)));

  if(tip == 0) {   /* size = 1 */
     *qbet = (int)image[atx][aty];
     *qalf = zeroalfa;
     *xd = 0;
     *yd = 0;
     *is = IDENTITY;
     return 0.0;
  }

  s0 = size * size;
  for(x=0;x < size;x++)
  for(y=0;y < size;y++) {
     pixel = (double)image[atx+x][aty+y];
     t0 +=  pixel;
     t2 +=  pixel * pixel;
     range[x][y] = pixel;
  }
  mean = t0 / s0;

  newclass(size,range,&isom,&clas);
  flips(size,range,flip_range,isom);
  stdd_class = std_class(size,flip_range);
  if(stdd_class > final_max_std) stdd_class = final_max_std;
  
  if (full_first_class) {
      start_first = 0;
      end_first   = 3;
  }
  else { 
      start_first = clas;
      end_first   = clas + 1;
  } 

  if (full_second_class) {
      start_second = 0;
      end_second   = 24; 
  }
  else { 
      start_second = stdd_class;
      end_second   = stdd_class + 1;
  } 

  for(fisher_first = start_first; fisher_first < end_first; fisher_first ++) 
  for(fisher_second=start_second; fisher_second < end_second; fisher_second++) { 
     pointer = class_std[tip][fisher_first][fisher_second];
     // printf("tip = %d, fisher_first=%d, fisher_second=%d\n",tip,fisher_first,fisher_second); 
     while(pointer != NULL) {
      // printf("inside while\n");
      // printf("isom : %d \t pointer->iso : %d\n",isom,pointer->iso);
      if(pointer->iso != NULL)
        isometry = mapping[isom][pointer->iso];
        comparisons ++ ;
        counter ++ ;
        s1 = pointer->sum;
        s2 = pointer->sum2;
        sd_d = pointer->var;
        t1 = 0.0;
        i = pointer->ptr_x >> 1;
        j = pointer->ptr_y >> 1;
   
        switch(isometry) { 
         case IDENTITY   :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[x][y];
                 break;
         case R_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[y][size-x-1];
                 break;
         case L_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-y-1][x];
                 break;
         case ROTATE180  :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-x-1][size-y-1];
                 break;
         case R_VERTICAL :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[x][size-y-1];
                 break;
         case R_HORIZONTAL: 
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-x-1][y];
                 break;
         case F_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[y][x];
                 break;
         case S_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-y-1][size-x-1];
                 break;
        }

        /* Compute the scalig factor */
        det = s2*s0 - s1*s1;  
        if(det == 0.0)
           alfa = 0.0;
        else
          alfa = (s0*t1 - s1*t0)  / det;
        if(alfa < 0.0 ) alfa = 0.0;

        /* Quantize the scalig factor */
        qalfa =  0.5 +(alfa )/( MAX_ALFA)* (1 << N_BITALFA);
        if(qalfa < 0)  qalfa = 0;
        if(qalfa >= (1 << N_BITALFA))  qalfa = (1 << N_BITALFA) - 1;

        /* Compute the scalig factor back from the quantized value*/
        alfa = (double) qalfa / (double)(1 << N_BITALFA) * ( MAX_ALFA) ;
        // printf("after deqauantization alfa = %f\n\n",alfa);

        /* Compute the offset */
        beta= (t0 - alfa*s1) / s0;
        if (alfa > 0.0)  beta += alfa * 255;

        /* Quantize the offset */
        qbeta = 0.5 + beta/ ((1.0+fabs(alfa))*255)*((1<< N_BITBETA)-1);
        if (qbeta< 0) qbeta = 0;
        if (qbeta >= 1 << N_BITBETA) qbeta = (1 << N_BITBETA)-1;
        // printf("%d\n",qbeta);
        /* Compute the offset back from the quantized value*/
        beta = (double)qbeta/(double)((1 << N_BITBETA)-1)*((1.0+fabs(alfa))*255);
        if (alfa > 0.0) beta  -= alfa * 255;

        /* Compute the rms distance */
        sum = t2 - 2*alfa*t1 -2*beta*t0 + alfa*alfa*s2 +
                              2*alfa*beta*s1 + s0*beta*beta;
        dist = sqrt(sum / s0);
         
        if(dist < min ) {
            min = dist;
            error_min = dist;
            alfa_min = alfa;
            beta_min = beta;
            *xd = pointer->ptr_x;
            *yd = pointer->ptr_y;
            *is = isometry;
            *qalf = qalfa;
            *qbet = qbeta;
        }
        pointer = pointer -> next;
        // printf("%f\nN",error_min);

     }
  }

// printf("%f \t %f \t %f\n",alfa_min,beta_min,error_min);
  // printf("Zero alfa count = %d\n",zero_alfa_count);
 return (min) ;
}

double BasicFIC_Coding(int atx,int aty,int size,int *xd,int *yd,int *is,
                                                       int* qalf,int *qbet)
{
  int x,y,qalfa ,qbeta,i,j;
  int tip,counter,isometry ;
  int isom,clas, var_class;
  int start_first, end_first, fisher_first;
  int start_second, end_second,fisher_second;
  double dist,alfa,beta, sd_r, sd_d;
  double min = 1000000000.0;
  double sum,s0;
  double mean,det;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  register double pixel;
  struct c *pointer;
  double alfa_min =0.0, beta_min=0.0, error_min = 0.0;
  tip = (int) rint((log((double) size)/log (2.0)));

  if(tip == 0) {   /* size = 1 */
     *qbet = (int)image[atx][aty];
     *qalf = zeroalfa;
     *xd = 0;
     *yd = 0;
     *is = IDENTITY;
     return 0.0;
  }

  s0 = size * size;
  for(x=0;x < size;x++)
  for(y=0;y < size;y++) {
     pixel = (double)image[atx+x][aty+y];
     t0 +=  pixel;
     t2 +=  pixel * pixel;
     range[x][y] = pixel;
  }
  mean = t0 / s0;

  newclass(size,range,&isom,&clas);
  flips(size,range,flip_range,isom);

  if (full_first_class) {
      start_first = 0;
      end_first   = 3;
  }
  else { 
      start_first = clas;
      end_first   = clas + 1;
  } 

  for(fisher_first = start_first; fisher_first < end_first; fisher_first ++){ 
     pointer = class_basicFIC[tip][fisher_first]; 
     while(pointer != NULL) {
        isometry = mapping[isom][pointer->iso];
        comparisons ++ ;
        counter ++ ;
        s1 = pointer->sum;
        s2 = pointer->sum2;
        sd_d = pointer->var;
        t1 = 0.0;

        i = pointer->ptr_x >> 1;
        j = pointer->ptr_y >> 1;
   
        switch(isometry) { 
         case IDENTITY   :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[x][y];
                 break;
         case R_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[y][size-x-1];
                 break;
         case L_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-y-1][x];
                 break;
         case ROTATE180  :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-x-1][size-y-1];
                 break;
         case R_VERTICAL :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[x][size-y-1];
                 break;
         case R_HORIZONTAL: 
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-x-1][y];
                 break;
         case F_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[y][x];
                 break;
         case S_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-y-1][size-x-1];
                 break;
        }

        /* Compute the scalig factor */
        det = s2*s0 - s1*s1;
        if(det == 0.0)
           alfa = 0.0;
        else
          alfa = (s0*t1 - s1*t0)  / det;
        if(alfa < 0.0 ) alfa = 0.0;
         // printf("before quantization alfa = %f\n",alfa);
        /* Quantize the scalig factor */
        // qalfa = 0.5 + (alfa )/( MAX_ALFA)* (1 << N_BITALFA);
        qalfa =  0.5 +(alfa )/( MAX_ALFA)* (1 << N_BITALFA);
        if(qalfa < 0)  qalfa = 0;
        if(qalfa >= (1 << N_BITALFA))  qalfa = (1 << N_BITALFA) - 1;

        /* Compute the scalig factor back from the quantized value*/
        alfa = (double) qalfa / (double)(1 << N_BITALFA) * ( MAX_ALFA) ;
//      
        if(alfa == 0)
          zero_alfa_count++;
        // printf("after deqauantization alfa = %f\n\n",alfa);

        /* Compute the offset */
        beta= (t0 - alfa*s1) / s0;
        if (alfa > 0.0)  beta += alfa * 255;

        /* Quantize the offset */
        qbeta = 0.5 + beta/ ((1.0+fabs(alfa))*255)*((1<< N_BITBETA)-1);
        if (qbeta< 0) qbeta = 0;
        if (qbeta >= 1 << N_BITBETA) qbeta = (1 << N_BITBETA)-1;
        // printf("%d\n",qbeta);
        /* Compute the offset back from the quantized value*/
        beta = (double)qbeta/(double)((1 << N_BITBETA)-1)*((1.0+fabs(alfa))*255);
        if (alfa > 0.0) beta  -= alfa * 255;

        /* Compute the rms distance */
        sum = t2 - 2*alfa*t1 -2*beta*t0 + alfa*alfa*s2 +
                              2*alfa*beta*s1 + s0*beta*beta;
        dist = sqrt(sum / s0);
         
        if(dist < min ) {
            min = dist;
            error_min = dist;
            alfa_min = alfa;
            beta_min = beta;
            *xd = pointer->ptr_x;
            *yd = pointer->ptr_y;
            *is = isometry;
            *qalf = qalfa;
            *qbet = qbeta;
        }
        pointer = pointer -> next;
        // printf("%f\n",error_min);

     }
  }

// printf("%f \t %f \t %f\n",alfa_min,beta_min,error_min);
  // printf("Zero alfa count = %d\n",zero_alfa_count);
 return (min) ;
}

double modified_FisherCoding_2(int atx,int aty,int size,int *xd,int *yd,int *is,
                                                       int* qalf,int *qbet)
{
  int x,y,qalfa ,qbeta,i,j;
  int tip,counter,isometry ;
  int isom,clas, var_class;
  int start_first, end_first, fisher_first;
  int start_second, end_second,fisher_second;
  double dist,alfa,beta, sd_r, sd_d;
  double min = 1000000000.0;
  double sum,s0;
  double mean,det;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  register double pixel;
  struct c *pointer;
  double alfa_min =0.0, beta_min=0.0, error_min = 0.0;
  tip = (int) rint((log((double) size)/log (2.0)));

  if(tip == 0) {   /* size = 1 */
     *qbet = (int)image[atx][aty];
     *qalf = zeroalfa;
     *xd = 0;
     *yd = 0;
     *is = IDENTITY;
     return 0.0;
  }

  s0 = size * size;
  for(x=0;x < size;x++)
  for(y=0;y < size;y++) {
     pixel = (double)image[atx+x][aty+y];
     t0 +=  pixel;
     t2 +=  pixel * pixel;
     range[x][y] = pixel;
  }
  mean = t0 / s0;
  sd_r = sqrt(variance_2(size,range,0,0));
  newclass(size,range,&isom,&clas);
  flips(size,range,flip_range,isom);
  var_class = variance_class(size,flip_range);

  if (full_first_class) {
      start_first = 0;
      end_first   = 3;
  }
  else { 
      start_first = clas;
      end_first   = clas + 1;
  } 

  if (full_second_class) {
      start_second = 0;
      end_second   = 24; 
  }
  else { 
      start_second = var_class;
      end_second   = var_class + 1;
  } 

  for(fisher_first = start_first; fisher_first < end_first; fisher_first ++) 
  for(fisher_second=start_second; fisher_second < end_second; fisher_second++) { 
     pointer = class_fisher[tip][fisher_first][fisher_second]; 
     while(pointer != NULL) {
        isometry = mapping[isom][pointer->iso];
        comparisons ++ ;
        counter ++ ;
        s1 = pointer->sum;
        s2 = pointer->sum2;
        t1 = 0.0;
        i = pointer->ptr_x >> 1;
        j = pointer->ptr_y >> 1;
        sd_d = sqrt(pointer->var);
        switch(isometry) { 
         case IDENTITY   :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[x][y];
                 break;
         case R_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[y][size-x-1];
                 break;
         case L_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-y-1][x];
                 break;
         case ROTATE180  :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-x-1][size-y-1];
                 break;
         case R_VERTICAL :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[x][size-y-1];
                 break;
         case R_HORIZONTAL: 
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-x-1][y];
                 break;
         case F_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[y][x];
                 break;
         case S_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-y-1][size-x-1];
                 break;
        }

        /* Compute the scalig factor */
        // det = s2*s0 - s1*s1;
        det = sd_d;
        if(det == 0.0)
           alfa = 0.0;
        else
          // alfa = (s0*t1 - s1*t0)  / det;
          alfa = sd_r /det;
        if(alfa < 0.0 ) alfa = 0.0;
         // printf("before quantization alfa = %f\n",alfa);
        /* Quantize the scalig factor */
        // qalfa = 0.5 + (alfa )/( MAX_ALFA)* (1 << N_BITALFA);
        qalfa =  0.5 +(alfa )/( MAX_ALFA)* (1 << N_BITALFA);
        if(qalfa < 0)  qalfa = 0;
        if(qalfa >= (1 << N_BITALFA))  qalfa = (1 << N_BITALFA) - 1;

        /* Compute the scalig factor back from the quantized value*/
        alfa = (double) qalfa / (double)(1 << N_BITALFA) * ( MAX_ALFA) ;
//      
        if(alfa == 0)
          zero_alfa_count++;
        // printf("after deqauantization alfa = %f\n\n",alfa);
       if((sd_r - (alfa * sd_d)) < (0.6 * sd_r)){    

          /* Compute the offset */
          beta= (t0 - alfa*s1) / s0;
          if (alfa > 0.0)  beta += alfa * 255;

          /* Quantize the offset */
          qbeta = 0.5 + beta/ ((1.0+fabs(alfa))*255)*((1<< N_BITBETA)-1);
          if (qbeta< 0) qbeta = 0;
          if (qbeta >= 1 << N_BITBETA) qbeta = (1 << N_BITBETA)-1;
          // printf("%d\n",qbeta);
          /* Compute the offset back from the quantized value*/
          beta = (double)qbeta/(double)((1 << N_BITBETA)-1)*((1.0+fabs(alfa))*255);
          if (alfa > 0.0) beta  -= alfa * 255;
   

          /* Compute the rms distance */
          sum = t2 - 2*alfa*t1 -2*beta*t0 + alfa*alfa*s2 +
                                2*alfa*beta*s1 + s0*beta*beta;
          dist = sqrt(sum / s0);
          // printf("alfa = %f \t beta = %f\n",alfa,beta);

          /* adaptive error based on range size */
          // sum = alfa * ( alfa*s2 - 2*t1 ) + t2 + qbeta * s0 * (qbeta - 2*beta);
          // dist  = sqrt(sum/s0);
          // printf("error = %f\n",sum);
          if(dist < min) {
              // printf("alfa = %ff\n",alfa);
              range_error = dist;
              min = dist;
              *xd = pointer->ptr_x;
              *yd = pointer->ptr_y;
              *is = isometry;
              *qalf = qalfa;
              *qbet = qbeta;
          }
        }
        pointer = pointer -> next;

     }
  }

 return (min) ;
}

double modified_FisherCoding_1(int atx,int aty,int size,int *xd,int *yd,int *is,
                                                       int* qalf,int *qbet)
{
  int x,y,qalfa ,qbeta,i,j;
  int tip,counter,isometry ;
  int isom,clas, var_class;
  int start_first, end_first, fisher_first;
  int start_second, end_second,fisher_second;
  double dist,alfa,beta, sd_r, sd_d;
  double min = 1000000000.0;
  double sum,s0;
  double mean,det;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  register double pixel;
  struct c *pointer;
  double alfa_min =0.0, beta_min=0.0, error_min = 0.0;
  tip = (int) rint((log((double) size)/log (2.0)));

  if(tip == 0) {   /* size = 1 */
     *qbet = (int)image[atx][aty];
     *qalf = zeroalfa;
     *xd = 0;
     *yd = 0;
     *is = IDENTITY;
     return 0.0;
  }

  s0 = size * size;
  for(x=0;x < size;x++)
  for(y=0;y < size;y++) {
     pixel = (double)image[atx+x][aty+y];
     t0 +=  pixel;
     t2 +=  pixel * pixel;
     range[x][y] = pixel;
  }
  mean = t0 / s0;
  sd_r = sqrt(variance_2(size,range,0,0));
  newclass(size,range,&isom,&clas);
  flips(size,range,flip_range,isom);
  var_class = variance_class(size,flip_range);

  if (full_first_class) {
      start_first = 0;
      end_first   = 3;
  }
  else { 
      start_first = clas;
      end_first   = clas + 1;
  } 

  if (full_second_class) {
      start_second = 0;
      end_second   = 24; 
  }
  else { 
      start_second = var_class;
      end_second   = var_class + 1;
  } 

  for(fisher_first = start_first; fisher_first < end_first; fisher_first ++) 
  for(fisher_second=start_second; fisher_second < end_second; fisher_second++) { 
     pointer = class_fisher[tip][fisher_first][fisher_second]; 
     while(pointer != NULL) {
        isometry = mapping[isom][pointer->iso];
        comparisons ++ ;
        counter ++ ;
        s1 = pointer->sum;
        s2 = pointer->sum2;
        t1 = 0.0;
        i = pointer->ptr_x >> 1;
        j = pointer->ptr_y >> 1;
        sd_d = sqrt(pointer->var);
        switch(isometry) { 
         case IDENTITY   :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[x][y];
                 break;
         case R_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[y][size-x-1];
                 break;
         case L_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-y-1][x];
                 break;
         case ROTATE180  :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-x-1][size-y-1];
                 break;
         case R_VERTICAL :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[x][size-y-1];
                 break;
         case R_HORIZONTAL: 
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-x-1][y];
                 break;
         case F_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[y][x];
                 break;
         case S_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-y-1][size-x-1];
                 break;
        }

        /* Compute the scalig factor */
        // det = s2*s0 - s1*s1;
        det = sd_d;
        if(det == 0.0)
           alfa = 0.0;
        else
          // alfa = (s0*t1 - s1*t0)  / det;
          alfa = sd_r /det;
        if(alfa < 0.0 ) alfa = 0.0;
         // printf("before quantization alfa = %f\n",alfa);
        /* Quantize the scalig factor */
        // qalfa = 0.5 + (alfa )/( MAX_ALFA)* (1 << N_BITALFA);
        qalfa =  0.5 +(alfa )/( MAX_ALFA)* (1 << N_BITALFA);
        if(qalfa < 0)  qalfa = 0;
        if(qalfa >= (1 << N_BITALFA))  qalfa = (1 << N_BITALFA) - 1;

        /* Compute the scalig factor back from the quantized value*/
        alfa = (double) qalfa / (double)(1 << N_BITALFA) * ( MAX_ALFA) ;
//      
        if(alfa == 0)
          zero_alfa_count++;
        // printf("after deqauantization alfa = %f\n\n",alfa);
       // if((sd_r - (alfa * sd_d)) < (0.6 * sd_r)){    

          /* Compute the offset */
          beta= (t0 - alfa*s1) / s0;
          if (alfa > 0.0)  beta += alfa * 255;

          /* Quantize the offset */
          qbeta = 0.5 + beta/ ((1.0+fabs(alfa))*255)*((1<< N_BITBETA)-1);
          if (qbeta< 0) qbeta = 0;
          if (qbeta >= 1 << N_BITBETA) qbeta = (1 << N_BITBETA)-1;
          // printf("%d\n",qbeta);
          /* Compute the offset back from the quantized value*/
          beta = (double)qbeta/(double)((1 << N_BITBETA)-1)*((1.0+fabs(alfa))*255);
          if (alfa > 0.0) beta  -= alfa * 255;
   

          /* Compute the rms distance */
          sum = t2 - 2*alfa*t1 -2*beta*t0 + alfa*alfa*s2 +
                                2*alfa*beta*s1 + s0*beta*beta;
          dist = sqrt(sum / s0);
          // printf("alfa = %f \t beta = %f\n",alfa,beta);

          /* adaptive error based on range size */
          // sum = alfa * ( alfa*s2 - 2*t1 ) + t2 + qbeta * s0 * (qbeta - 2*beta);
          // dist  = sqrt(sum/s0);
          // printf("error = %f\n",sum);
          if(dist < min) {
              // printf("alfa = %ff\n",alfa);
              range_error = dist;
              min = dist;
              *xd = pointer->ptr_x;
              *yd = pointer->ptr_y;
              *is = isometry;
              *qalf = qalfa;
              *qbet = qbeta;
          }
        //}
        pointer = pointer -> next;

     }
  }

 return (min) ;
}

double modified_FisherCoding(int atx,int aty,int size,int *xd,int *yd,int *is,
                                                       int* qalf,int *qbet)
{
  int x,y,qalfa ,qbeta,i,j;
  int tip,counter,isometry ;
  int isom,clas, var_class;
  int start_first, end_first, fisher_first;
  int start_second, end_second,fisher_second;
  double dist,alfa,beta, sd_r, sd_d;
  double min = 1000000000.0;
  double sum,s0;
  double mean,det;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  register double pixel;
  struct c *pointer;
  double alfa_min =0.0, beta_min=0.0, error_min = 0.0;
  tip = (int) rint((log((double) size)/log (2.0)));

  if(tip == 0) {   /* size = 1 */
     *qbet = (int)image[atx][aty];
     *qalf = zeroalfa;
     *xd = 0;
     *yd = 0;
     *is = IDENTITY;
     return 0.0;
  }

  s0 = size * size;
  for(x=0;x < size;x++)
  for(y=0;y < size;y++) {
     pixel = (double)image[atx+x][aty+y];
     t0 +=  pixel;
     t2 +=  pixel * pixel;
     range[x][y] = pixel;
  }
  mean = t0 / s0;
  sd_r = sqrt(variance_2(size,range,0,0));
  newclass(size,range,&isom,&clas);
  flips(size,range,flip_range,isom);
  var_class = variance_class(size,flip_range);

  if (full_first_class) {
      start_first = 0;
      end_first   = 3;
  }
  else { 
      start_first = clas;
      end_first   = clas + 1;
  } 

  if (full_second_class) {
      start_second = 0;
      end_second   = 24; 
  }
  else { 
      start_second = var_class;
      end_second   = var_class + 1;
  } 

  for(fisher_first = start_first; fisher_first < end_first; fisher_first ++) 
  for(fisher_second=start_second; fisher_second < end_second; fisher_second++) { 
     pointer = class_fisher[tip][fisher_first][fisher_second]; 
     while(pointer != NULL) {
        isometry = mapping[isom][pointer->iso];
        comparisons ++ ;
        counter ++ ;
        s1 = pointer->sum;
        s2 = pointer->sum2;
        t1 = 0.0;
        i = pointer->ptr_x >> 1;
        j = pointer->ptr_y >> 1;
        sd_d = sqrt(pointer->var);
        switch(isometry) { 
         case IDENTITY   :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[x][y];
                 break;
         case R_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[y][size-x-1];
                 break;
         case L_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-y-1][x];
                 break;
         case ROTATE180  :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-x-1][size-y-1];
                 break;
         case R_VERTICAL :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[x][size-y-1];
                 break;
         case R_HORIZONTAL: 
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-x-1][y];
                 break;
         case F_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[y][x];
                 break;
         case S_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-y-1][size-x-1];
                 break;
        }

        /* Compute the scalig factor */
        // det = s2*s0 - s1*s1;
        det = sd_d;
        if(det == 0.0)
           alfa = 0.0;
        else
          // alfa = (s0*t1 - s1*t0)  / det;
          alfa = sd_r /det;
        if(alfa < 0.0 ) alfa = 0.0;
         // printf("before quantization alfa = %f\n",alfa);
        /* Quantize the scalig factor */
        // qalfa = 0.5 + (alfa )/( MAX_ALFA)* (1 << N_BITALFA);
        qalfa =  0.5 +(alfa )/( MAX_ALFA)* (1 << N_BITALFA);
        if(qalfa < 0)  qalfa = 0;
        if(qalfa >= (1 << N_BITALFA))  qalfa = (1 << N_BITALFA) - 1;

        /* Compute the scalig factor back from the quantized value*/
        alfa = (double) qalfa / (double)(1 << N_BITALFA) * ( MAX_ALFA) ;
//      
        if(alfa == 0)
          zero_alfa_count++;
        // printf("after deqauantization alfa = %f\n\n",alfa);
       // if((sd_r - (alfa * sd_d)) < (0.6 * sd_r)){    

          /* Compute the offset */
          beta= (t0 - alfa*s1) / s0;
          if (alfa > 0.0)  beta += alfa * 255;

          /* Quantize the offset */
          qbeta = 0.5 + beta/ ((1.0+fabs(alfa))*255)*((1<< N_BITBETA)-1);
          if (qbeta< 0) qbeta = 0;
          if (qbeta >= 1 << N_BITBETA) qbeta = (1 << N_BITBETA)-1;
          // printf("%d\n",qbeta);
          /* Compute the offset back from the quantized value*/
          beta = (double)qbeta/(double)((1 << N_BITBETA)-1)*((1.0+fabs(alfa))*255);
          if (alfa > 0.0) beta  -= alfa * 255;
   

          /* Compute the rms distance */
          sum = t2 - 2*alfa*t1 -2*beta*t0 + alfa*alfa*s2 +
                                2*alfa*beta*s1 + s0*beta*beta;
          dist = sqrt(sum / s0);
          // printf("alfa = %f \t beta = %f\n",alfa,beta);

          /* adaptive error based on range size */
          // sum = alfa * ( alfa*s2 - 2*t1 ) + t2 + qbeta * s0 * (qbeta - 2*beta);
          // dist  = sqrt(sum/s0);
          // printf("error = %f\n",sum);
          if(dist < min) {
              // printf("alfa = %ff\n",alfa);
              range_error = dist;
              min = dist;
              *xd = pointer->ptr_x;
              *yd = pointer->ptr_y;
              *is = isometry;
              *qalf = qalfa;
              *qbet = qbeta;
          }
        //}
        pointer = pointer -> next;

     }
  }

 return (min) ;
}

double FisherCoding(int atx,int aty,int size,int *xd,int *yd,int *is,
                                                       int* qalf,int *qbet)
{
  int x,y,qalfa ,qbeta,i,j;
  int tip,counter,isometry ;
  int isom,clas, var_class;
  int start_first, end_first, fisher_first;
  int start_second, end_second,fisher_second;
  double dist,alfa,beta, sd_r, sd_d;
  double min = 1000000000.0;
  double sum,s0;
  double mean,det;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  register double pixel;
  struct c *pointer;
  double alfa_min =0.0, beta_min=0.0, error_min = 0.0;
  tip = (int) rint((log((double) size)/log (2.0)));

  if(tip == 0) {   /* size = 1 */
     *qbet = (int)image[atx][aty];
     *qalf = zeroalfa;
     *xd = 0;
     *yd = 0;
     *is = IDENTITY;
     return 0.0;
  }

  s0 = size * size;
  for(x=0;x < size;x++)
  for(y=0;y < size;y++) {
     pixel = (double)image[atx+x][aty+y];
     t0 +=  pixel;
     t2 +=  pixel * pixel;
     range[x][y] = pixel;
  }
  mean = t0 / s0;

  newclass(size,range,&isom,&clas);
  flips(size,range,flip_range,isom);
  var_class = variance_class(size,flip_range);

  if (full_first_class) {
      start_first = 0;
      end_first   = 3;
  }
  else { 
      start_first = clas;
      end_first   = clas + 1;
  } 

  if (full_second_class) {
      start_second = 0;
      end_second   = 24; 
  }
  else { 
      start_second = var_class;
      end_second   = var_class + 1;
  } 

  for(fisher_first = start_first; fisher_first < end_first; fisher_first ++) 
  for(fisher_second=start_second; fisher_second < end_second; fisher_second++) { 
     pointer = class_fisher[tip][fisher_first][fisher_second]; 
     while(pointer != NULL) {
        isometry = mapping[isom][pointer->iso];
        comparisons ++ ;
        counter ++ ;
        s1 = pointer->sum;
        s2 = pointer->sum2;
        t1 = 0.0;
        i = pointer->ptr_x >> 1;
        j = pointer->ptr_y >> 1;
   
        switch(isometry) { 
         case IDENTITY   :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[x][y];
                 break;
         case R_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[y][size-x-1];
                 break;
         case L_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-y-1][x];
                 break;
         case ROTATE180  :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-x-1][size-y-1];
                 break;
         case R_VERTICAL :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[x][size-y-1];
                 break;
         case R_HORIZONTAL: 
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-x-1][y];
                 break;
         case F_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[y][x];
                 break;
         case S_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j]*range[size-y-1][size-x-1];
                 break;
        }

        /* Compute the scalig factor */
        det = s2*s0 - s1*s1;
        // det = sd_d;
        if(det == 0.0)
           alfa = 0.0;
        else
          alfa = (s0*t1 - s1*t0)  / det;
          // alfa = sd_r /det;
        if(alfa < 0.0 ) alfa = 0.0;
         // printf("before quantization alfa = %f\n",alfa);
        /* Quantize the scalig factor */
        // qalfa = 0.5 + (alfa )/( MAX_ALFA)* (1 << N_BITALFA);
        qalfa =  0.5 +(alfa )/( MAX_ALFA)* (1 << N_BITALFA);
        if(qalfa < 0)  qalfa = 0;
        if(qalfa >= (1 << N_BITALFA))  qalfa = (1 << N_BITALFA) - 1;

        /* Compute the scalig factor back from the quantized value*/
        alfa = (double) qalfa / (double)(1 << N_BITALFA) * ( MAX_ALFA) ;
//      
        if(alfa == 0)
          zero_alfa_count++;
        // printf("after deqauantization alfa = %f\n\n",alfa);

        /* Compute the offset */
        beta= (t0 - alfa*s1) / s0;
        if (alfa > 0.0)  beta += alfa * 255;

        /* Quantize the offset */
        qbeta = 0.5 + beta/ ((1.0+fabs(alfa))*255)*((1<< N_BITBETA)-1);
        if (qbeta< 0) qbeta = 0;
        if (qbeta >= 1 << N_BITBETA) qbeta = (1 << N_BITBETA)-1;
        // printf("%d\n",qbeta);
        /* Compute the offset back from the quantized value*/
        beta = (double)qbeta/(double)((1 << N_BITBETA)-1)*((1.0+fabs(alfa))*255);
        if (alfa > 0.0) beta  -= alfa * 255;
 

        /* Compute the rms distance */
        sum = t2 - 2*alfa*t1 -2*beta*t0 + alfa*alfa*s2 +
                              2*alfa*beta*s1 + s0*beta*beta;
        dist = sqrt(sum / s0);
        // printf("alfa = %f \t beta = %f\n",alfa,beta);

        /* adaptive error based on range size */
        // sum = alfa * ( alfa*s2 - 2*t1 ) + t2 + qbeta * s0 * (qbeta - 2*beta);
        // dist  = sqrt(sum/s0);
        // printf("error = %f\n",sum);
        if(dist < min/*max_adaptive_error * size * size*/ ) {
            // printf("alfa = %ff\n",alfa);
            range_error = dist;
            min = dist;
            *xd = pointer->ptr_x;
            *yd = pointer->ptr_y;
            *is = isometry;
            *qalf = qalfa;
            *qbet = qbeta;
        }
        pointer = pointer -> next;

     }
  }

 return (min) ;
}

int getL1class(int x_size, int y_size, int old_r_class)
{
  if(adaptive_fisher_class[y_size][x_size][old_r_class] != NULL){
    return old_r_class;
  }
  
  int new_r_class = old_r_class;
  int L1_length = 24;
  
  for(int fwd = 0 , bwd = L1_length - 1;fwd <= L1_length - 1, bwd >= 0;){
     

     if(adaptive_fisher_class[y_size][x_size][fwd] != NULL ){
      return fwd;
     }
     else if(adaptive_fisher_class[y_size][x_size][bwd] != NULL ){
      return bwd;
     }
     fwd++;
     bwd--;

  }

 
  return new_r_class;
}

int quantize(double value, double max, int imax)
{
  int ival = (int)(value*imax/max);

  if (ival < 0)
  {
    return 0;
  }

  if (ival > imax)
  {
    return imax;
  }

  return ival;
}



double dequantize(double value, double max, double imax)
{
  return ((double)(value)*(max)/(double)imax);
}


double adaptiveFisherCoding(int atx,int aty,int x_size,int y_size,int *xd,int *yd,int *is,
                                                       int* qalf,int *qbet)
{
  int x,y,i,j;
  int tip_x,tip_y,counter,isometry ;
  int isom,clas, var_class;
  int start_first, end_first, fisher_first;
  int start_second, end_second,fisher_second;
  double dist,alfa,beta, sd_r, sd_d;
  double min = 1000000000.0;
  double sum,s0;
  double mean,det;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  register double pixel;
  struct c *pointer;
  double alfa_min =0.0, beta_min=0.0, error_min = 0.0;
  tip_y =  y_size;
  tip_x =   x_size;
  int qalfa, qbeta;
     // printf("x_size = %d\t y_size = %d\n",x_size,y_size);

  // if(x_size == 1  && y_size == 1) {   /* size = 1 */
  //    // printf("x_size = %d\t y_size = %d\n",x_size,y_size);
  //    *qbet = (int)image[atx][aty];
  //    *qalf = zeroalfa;
  //    *xd = 0;
  //    *yd = 0;
  //    *is = IDENTITY;
  //    return 0.0;
  // }
  // printf("xsize= %d \t ysize = %d\n",x_size,y_size );
 
  s0 = y_size * x_size;
  for(int y=0;y < y_size;y++){
    for(int x=0;x < x_size;x++) {
       pixel = (double)image[aty+y][atx+x];
       t0 +=  pixel;
       t2 +=  pixel * pixel;
       // printf("x= %d \t y = %d\n",x,y );
       range[y][x] = pixel;
    }
  }
  // sd_r = sqrt(variance_3(x_size,y_size, range, 0,0));
  mean = t0 / s0;
      
  adaptiveNewclass_2(x_size,y_size,range,&isom,&clas); /* brings the domain in the canonical   */
  // printf("old Range class  = %d\n",clas);
  // flips(size,range,flip_range,isom);
  //var_class = variance_class(size,flip_range);
  // clas = getL1class(2*x_size,2*y_size,clas);
  // printf("new range class = %d\n",clas);
  // if(adaptive_fisher_class[y_size][x_size][clas] == NULL)
    // printf("adaptive_fisher_class = NULL");

  if (full_first_class) {
      start_first = 0;
      end_first   = 24;
  }
  else { 
      start_first = clas;
      end_first   = clas +1 ;
  } 

  // if (full_second_class) {
  //     start_second = 0;
  //     end_second   = 24; 
  // }
  // else { 
  //     start_second = var_class;
  //     end_second   = var_class + 1;
  // } 

  for(fisher_first = start_first; fisher_first < end_first; fisher_first++){ 
   // for(fisher_second=start_second; fisher_second < end_second; fisher_second++){ 
     pointer = adaptive_fisher_class[tip_y][tip_x][fisher_first]; 
     while(pointer != NULL) {
        // isometry = mapping[isom][pointer->iso];
        comparisons ++ ;
        counter ++ ;
        s1 = pointer->sum;
        s2 = pointer->sum2;
        sd_d = pointer->var;
        t1 = 0.0;
        i = pointer->ptr_x >> 1;
        j = pointer->ptr_y >> 1;
        isometry = 0;
        sum = 0.0;
        alfa = beta = 0.0;
        // printf("contract[214][216] = %u\n",contract[200][256]);
        switch(isometry) { 
         case IDENTITY   :  
                 for(x=0;x< y_size;x++)
                 for(y=0;y< x_size;y++)
                     t1 += contract[x+i][y+j]*range[x][y];
                 break;
         // case R_ROTATE90 :  
         //         for(x=0;x< y_size;x++)
         //         for(y=0;y< x_size;y++)
         //             t1 += contract[x+i][y+j]*range[y][y_size-x-1];
         //         break;
         // case L_ROTATE90 :  
         //         for(x=0;x< y_size;x++)
         //         for(y=0;y< x_size;y++)
         //             t1 += contract[x+i][y+j]*range[x_size-y-1][x];
         //         break;
         // case ROTATE180  :  
         //         for(x=0;x< y_size;x++)
         //         for(y=0;y< x_size;y++)
         //             t1 += contract[x+i][y+j]*range[y_size-x-1][x_size-y-1];
         //         break;
         // case R_VERTICAL :  
         //         for(x=0;x< y_size;x++)
         //         for(y=0;y< x_size;y++)
         //             t1 += contract[x+i][y+j]*range[x][x_size-y-1];
         //         break;
         // case R_HORIZONTAL: 
         //         for(x=0;x< y_size;x++)
         //         for(y=0;y< x_size;y++)
         //             t1 += contract[x+i][y+j]*range[y_size-x-1][y];
         //         break;
         // case F_DIAGONAL:   
         //         for(x=0;x< y_size;x++)
         //         for(y=0;y< x_size;y++)
         //             t1 += contract[x+i][y+j]*range[y][x];
         //         break;
         // case S_DIAGONAL:   
         //         for(x=0;x< y_size;x++)
         //         for(y=0;y< x_size;y++)
         //             t1 += contract[x+i][y+j]*range[x_size-y-1][y_size-x-1];
         //         break;
        }
        // printf("X size = %d \t Y size = %d\n",x_size,y_size);
        // printf("R sum = %f\n",t0);
        // printf("R sum2 = %f\n",t2);
        // printf("D sum = %f\n",s1);
        // printf("D sum2 = %f\n",s2);

        double original_beta;
        /* Compute the scalig factor */
        det = s2*s0 - s1*s1;
        if(det == 0.0)
           alfa = 0.0; 
        else
          alfa = (s0*t1 - s1*t0)  / det;
        if(alfa < 0.0 ) alfa = 0.0;
        alfa = 0.0;
        /* Quantize the scalig factor */
        MAX_ALFA = 1.0;
        N_BITALFA = 5;
        N_BITBETA = 7;
        qalfa = quantize(alfa, MAX_ALFA, (int)(1 << N_BITALFA) - 1 );

        alfa = dequantize(qalfa,MAX_ALFA,(int)(1 << N_BITALFA) - 1 );
        /* compute offset */
        beta = (t0 - alfa * s1) / s0;    
        original_beta = beta;
        double max_scaled = (double)(alfa * 255.0);
        qbeta = quantize(beta + max_scaled, max_scaled + 255, (int)(1 << N_BITBETA));
        // qbeta = 2 * qbeta;
        beta = dequantize(qbeta,max_scaled + 255, (int)(1 << N_BITBETA));
        
        /* Compute the rms distance */
        // sum = t2 - 2*alfa*t1 -2*beta*t0 + alfa*alfa*s2 +
        //                       2*alfa*beta*s1 + s0*beta*beta;
        // dist = sqrt(sum / s0);
         // printf("alfa = %f\n",alfa);
         // printf("beta = %f\n",beta);

        sum = alfa * (alfa * s2 - 2*t1) + t2 + 
                     (double)beta * s0 * ((double)beta - 2.0*original_beta);  
         // printf("error2 = %f\n",sum);
        
        // printf("alfa = %f \t beta = %f\n",alfa,beta);

        /* adaptive error based on range size */
        // sum = alfa * ( alfa*s2 - 2*t1 ) + t2 + qbeta * s0 * (qbeta - 2*beta);
        // dist  = sqrt(sum/s0);
        // printf("error = %f\n",sum);
        if(sum < min ) {
            // printf("error2 = %f\n",sum);
            min = sum;
            *xd = pointer->ptr_x;
            *yd = pointer->ptr_y;
            *is = isometry;
            *qalf = qalfa;
            *qbet = qbeta;
        }
        pointer = pointer -> next;

     }
    //}
  }
   // printf("%f \t %f \t %f\n",alfa_min,beta_min,error_min);
        // printf("min=%f\n",min);

   // printf("Zero alfa count = %d\n",zero_alfa_count);
 return (min) ;
}


double HurtgenCoding(int atx,int aty,int size,int *xd,int *yd,int *is,
                                                       int* qalf,int *qbet)
{
  int x,y,qalfa ,qbeta,i,j;
  int tip,start_first, end_first;
  int start_second, end_second;
  int counter,isom;
  int m_class_1, m_class_2;
  int hur_first,hur_second;
  double dist,alfa,beta;
  double min = 1000000000.0;
  double sum,s0,mean,det;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  register double pixel;
  struct c *pointer;


  tip = (int)rint((log((double) size)/log (2.0)));

  if(tip == 0) {   /* size = 1 */
      *qbet = (int)image[atx][aty];
      *qalf = zeroalfa;
      *xd = 0;
      *yd = 0;
      *is = IDENTITY;
      return 0.0;
  }

  s0 = size * size;
  for(x=0;x < size;x++)
  for(y=0;y < size;y++) {
     pixel = (double)image[atx+x][aty+y];
     t0 +=  pixel;
     t2 +=  pixel * pixel;
     range[x][y] = pixel;
  }

  mean = t0 / s0;

  for(isom=IDENTITY;isom <= S_DIAGONAL;isom++) {

     switch(isom) {     
        case L_ROTATE90  : flips(size,range,flip_range,R_ROTATE90);
                           break;
        case R_ROTATE90  : flips(size,range,flip_range,L_ROTATE90);
                           break; 
        default          : flips(size,range,flip_range,isom);
                           break;
     }
     m_class_1 = hurtgen_class(size,flip_range);
     m_class_2 = variance_class(size,flip_range);

     if (full_first_class) {
        start_first = 0;
        end_first   = 15;   /* Class 15 is always empty */
     }
     else { 
        start_first = m_class_1;
        end_first   = m_class_1 + 1;
     } 

     if (full_second_class) {
        start_second = 0;
        end_second   = 24; 
     }
     else { 
        start_second = m_class_2;
        end_second   = m_class_2 + 1;
     } 

     for(hur_first = start_first; hur_first < end_first; hur_first ++) 
     for(hur_second = start_second; hur_second < end_second; hur_second++) { 
        pointer = class_hurtgen[tip][hur_first][hur_second];
        while(pointer != NULL) {
           comparisons ++ ;
           counter ++ ;
           s1 = pointer->sum;
           s2 = pointer->sum2;
           t1 = 0.0;
           i = pointer->ptr_x >> 1;
           j = pointer->ptr_y >> 1;
   
           switch(isom) {
            case IDENTITY   :  
                    for(x=0;x< size;x++)
                    for(y=0;y< size;y++)
                       t1 += contract[x+i][y+j] * range[x][y];
                    break;
            case R_ROTATE90 :  
                    for(x=0;x< size;x++)
                    for(y=0;y< size;y++)
                       t1 += contract[x+i][y+j] * range[y][size-x-1];
                    break;
            case L_ROTATE90 :  
		    for(x=0;x< size;x++)
                    for(y=0;y< size;y++)
                       t1 += contract[x+i][y+j] * range[size-y-1][x];
                    break;
            case ROTATE180  :  
                    for(x=0;x< size;x++)
                    for(y=0;y< size;y++)
                       t1 += contract[x+i][y+j] * range[size-x-1][size-y-1];
                    break;
            case R_VERTICAL :  
                    for(x=0;x< size;x++)
                    for(y=0;y< size;y++)
                       t1 += contract[x+i][y+j] * range[x][size-y-1];
                    break;
            case R_HORIZONTAL: 
                    for(x=0;x< size;x++)
                    for(y=0;y< size;y++)
                       t1 += contract[x+i][y+j] * range[size-x-1][y];
                    break;
            case F_DIAGONAL:   
                    for(x=0;x< size;x++)
                    for(y=0;y< size;y++)
                       t1 += contract[x+i][y+j] * range[y][x];
                    break;
            case S_DIAGONAL:   
                    for(x=0;x< size;x++)
                    for(y=0;y< size;y++)
                        t1 += contract[x+i][y+j] * range[size-y-1][size-x-1];
                    break;
           }

           /* Compute the scalig factor */
           det = s2*s0 - s1*s1;
           if(det == 0.0)
              alfa = 0.0;
           else
             alfa = (s0*t1 - s1*t0) / det;
           if(alfa < 0.0 ) alfa = 0.0;

           /* Quantize the scalig factor */
           qalfa = 0.5 + (alfa )/( MAX_ALFA)* (1 << N_BITALFA);
           if(qalfa < 0)  qalfa = 0;
           if(qalfa >= (1 << N_BITALFA))  qalfa = (1 << N_BITALFA) - 1;

           /* Compute the scalig factor back from the quantized value*/
           alfa = (double) qalfa / (double)(1 << N_BITALFA) * ( MAX_ALFA) ;

           /* Compute the offset */
           beta= (t0 - alfa*s1) / s0;
           if (alfa > 0.0)  beta += alfa * 255;

           /* Quantize the offset */
           qbeta = 0.5 + beta/ ((1.0+fabs(alfa))*255)*((1 << N_BITBETA)-1);
           if (qbeta< 0) qbeta = 0;
           if (qbeta >= 1 << N_BITBETA) qbeta = (1 << N_BITBETA)-1;

           /* Compute the offset back from the quantized value*/
           beta = (double) qbeta/(double)((1 << N_BITBETA)-1)*
                                             ((1.0+fabs(alfa)) * 255);
           if (alfa > 0.0) beta  -= alfa * 255;

           /* Compute the rms distance */
           sum = t2 - 2*alfa*t1 -2*beta*t0 + alfa*alfa*s2 +
                              2*alfa*beta*s1 + s0*beta*beta;
           dist = sqrt(sum / s0);

           if(dist < min ) {
              min = dist;
              *xd = pointer->ptr_x;
              *yd = pointer->ptr_y;
              *is = isom;
              *qalf = qalfa;
              *qbet = qbeta;
           }
           pointer = pointer -> next;
        }
     }
  }
  return (min) ;
}


double Mc_SaupeCoding(int atx,int aty,int size,int *xd,int *yd,int *is,
                                                         int* qalf,int *qbet)
{
  int x,y,qalfa,qbeta,i,j;
  int counter,found = 0;
  int tip,ii,isom,a,b,r,t;
  int aa,bb,cx,cy;
  double dist,alfa,beta;
  double min = 1000000000.0;
  double sum,s0,mean,det;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  double r_vector_mc[16];
  static float r_vector_saupe[4096];
  static int nlist[MAX_NEIGHBOURS]; 
  register double pixel;


  tip = (int)rint((log((double) size)/log (2.0)));

  if(tip == 0) {   /* size = 1 */
     *qbet = (int)image[atx][aty];
     *qalf = zeroalfa;
     *xd = 0;
     *yd = 0;
     *is = IDENTITY;
     return 0.0;
  }

  s0 = size * size;
  for(x=0;x < size;x++)
  for(y=0;y < size;y++) {
     pixel = (double)image[atx+x][aty+y];
     t0 +=  pixel;
     t2 +=  pixel * pixel;
     range[x][y] = pixel;
  }

  mean = t0 / s0;

  for(isom=IDENTITY;isom <= S_DIAGONAL;isom++) {
     switch(isom) {     
         case L_ROTATE90  : flips(size,range,flip_range,R_ROTATE90);
                            flips(size,range,range_tmp,R_ROTATE90);
                            break;
         case R_ROTATE90  : flips(size,range,flip_range,L_ROTATE90);
                            flips(size,range,range_tmp,L_ROTATE90);
                            break; 
         default          : flips(size,range,flip_range,isom);
                            flips(size,range,range_tmp,isom);
                            break;
     }

     ComputeSaupeVectors(range_tmp,size,tip,r_vector_saupe);
     ComputeMcVectors(flip_range,range_tmp,size,tip,r_vector_mc);

     cx = (int) (r_vector_mc[0] / TWOPI * n_p_class ); 
     cy = (int) (r_vector_mc[1] / TWOPI * n_p_class ); 

     r = 0;
     t = 1;

     while(1) {  /* To make sure there is at least one kd-tree searched  */
       for(a=(cx+n_p_class-r)%n_p_class, aa=0;aa<t;a=(a+1)%n_p_class, aa++) 
       for(b=(cy+n_p_class-r)%n_p_class, bb=0;bb<t;b=(b+1)%n_p_class, bb++) 
         if( class_polar[tip][a][b] != NULL)
            goto out_loops;

       r += 1;
       t += 1;
     }   
      
  out_loops:

     for(a=(cx+n_p_class-r)%n_p_class, aa=0;aa<t;a=(a+1)%n_p_class, aa++) 
     for(b=(cy+n_p_class-r)%n_p_class, bb=0;bb<t;b=(b+1)%n_p_class, bb++) {

        found = 0;
        if(class_polar_saupe[tip][a][b] != NULL)
            found = kdtree_search(r_vector_saupe, f_vect[tip][a][b], 
                       feat_vect_dim[tip], class_polar_saupe[tip][a][b], eps, 
                                                                    matches,nlist); 

        for(ii=0;ii<found;ii++) {
           comparisons ++ ;
           counter ++ ;
           s1 = c_book[tip][a][b][nlist[ii]].sum;
           s2 = c_book[tip][a][b][nlist[ii]].sum2;
           t1 = 0.0;
           i = c_book[tip][a][b][nlist[ii]].ptr_x >> 1;
           j = c_book[tip][a][b][nlist[ii]].ptr_y >> 1;
   
           switch(isom) {
            case IDENTITY   :  
                    for(x=0;x< size;x++)
                    for(y=0;y< size;y++)
                        t1 += contract[x+i][y+j] * range[x][y];
                    break;
            case R_ROTATE90 :  
                    for(x=0;x< size;x++)
                    for(y=0;y< size;y++)
                        t1 += contract[x+i][y+j] * range[y][size-x-1];
                    break;
            case L_ROTATE90 :  
                    for(x=0;x< size;x++)
                    for(y=0;y< size;y++)
                        t1 += contract[x+i][y+j] * range[size-y-1][x];
                    break;
            case ROTATE180  :  
                    for(x=0;x< size;x++)
                    for(y=0;y< size;y++)
                        t1 += contract[x+i][y+j] * range[size-x-1][size-y-1];
                    break;
            case R_VERTICAL :  
                    for(x=0;x< size;x++)
                    for(y=0;y< size;y++)
                        t1 += contract[x+i][y+j] * range[x][size-y-1];
                    break;
            case R_HORIZONTAL: 
                    for(x=0;x< size;x++)
                    for(y=0;y< size;y++)
                        t1 += contract[x+i][y+j] * range[size-x-1][y];
                    break;
            case F_DIAGONAL:   
                    for(x=0;x< size;x++)
                    for(y=0;y< size;y++)
                        t1 += contract[x+i][y+j] * range[y][x];
                    break;
            case S_DIAGONAL:   
                    for(x=0;x< size;x++)
                    for(y=0;y< size;y++)
                        t1 += contract[x+i][y+j] * range[size-y-1][size-x-1];
                    break;
           }

           /* Compute the scalig factor */
           det = s2*s0 - s1*s1;
           if(det == 0.0)
              alfa = 0.0;
           else
             alfa = (s0*t1 - s1*t0) / det;
           if(alfa < 0.0 ) alfa = 0.0;

           /* Quantize the scalig factor */
           qalfa = 0.5 + (alfa )/( MAX_ALFA)* (1 << N_BITALFA);
           if(qalfa < 0)  qalfa = 0;
           if(qalfa >= (1 << N_BITALFA))  qalfa = (1 << N_BITALFA) - 1;

           /* Compute the scalig factor back from the quantized value*/
           alfa = (double) qalfa / (double)(1 << N_BITALFA) * ( MAX_ALFA) ;

           /* Compute the offset */
           beta= (t0 - alfa*s1) / s0;
           if (alfa > 0.0)  beta += alfa * 255;

           /* Quantize the offset */
           qbeta = 0.5 + beta/ ((1.0+fabs(alfa))*255)*((1 << N_BITBETA)-1);
           if (qbeta< 0) qbeta = 0;
           if (qbeta >= 1 << N_BITBETA) qbeta = (1 << N_BITBETA)-1;

           /* Compute the offset back from the quantized value*/
           beta = (double) qbeta/(double)((1 << N_BITBETA)-1)*
                                                 ((1.0+fabs(alfa)) * 255);
           if (alfa > 0.0) beta  -= alfa * 255;

           /* Compute the rms distance */
           sum = t2 - 2*alfa*t1 -2*beta*t0 + alfa*alfa*s2 +
                              2*alfa*beta*s1 + s0*beta*beta;
           dist = sqrt(sum / s0);

           if(dist < min ) {
              min = dist;
              *xd = c_book[tip][a][b][nlist[ii]].ptr_x;
              *yd = c_book[tip][a][b][nlist[ii]].ptr_y; 
              *is = isom;
              *qalf = qalfa;
              *qbet = qbeta;
           }
        }
     }
  }
  return (min) ;
}


double MassCenterCoding(int atx,int aty,int size,int *xd,int *yd,
                                                int *is, int* qalf,int *qbet)
{
  int x,y,qalfa ,qbeta,i,j;
  int tip,counter,isom;
  int a,b,r,t,aa,bb,cx,cy;
  double dist,alfa,beta;
  double min = 1000000000.0;
  double sum,s0,mean,det;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  double r_vector[16]; 
  register double pixel;
  struct c *pointer;


  tip = (int)rint((log((double) size)/log (2.0)));

  if(tip == 0) {   /* size = 1 */
      *qbet = (int)image[atx][aty];
      *qalf = zeroalfa;
      *xd = 0;
      *yd = 0;
      *is = IDENTITY;
      return 0.0;
  }


  s0 = size * size;
  for(x=0;x < size;x++)
  for(y=0;y < size;y++) {
     pixel = (double)image[atx+x][aty+y];
     t0 +=  pixel;
     t2 +=  pixel * pixel;
     range[x][y] = pixel;
  }

  mean = t0 / s0;

  for(isom=IDENTITY;isom <= S_DIAGONAL;isom++) {
     switch(isom) {     
         case L_ROTATE90  : flips(size,range,flip_range,R_ROTATE90);
                            break;
         case R_ROTATE90  : flips(size,range,flip_range,L_ROTATE90);
                            break; 
         default          : flips(size,range,flip_range,isom);
                            break;
      }

      ComputeMcVectors(flip_range,range_tmp,size,tip,r_vector);

      cx = (int) (r_vector[0] / TWOPI * n_p_class ); 
      cy = (int) (r_vector[1] / TWOPI * n_p_class ); 

      r = 1;
      t = 3;

      while(1) { /* To make sure there is at least one list searched */
        for(a=(cx+n_p_class-r)%n_p_class, aa=0;aa<t;a=(a+1)%n_p_class, aa++) 
        for(b=(cy+n_p_class-r)%n_p_class, bb=0;bb<t;b=(b+1)%n_p_class, bb++) 
          if( class_polar[tip][a][b] != NULL)
             goto out_loops;

        r += 1;
        t += 1;
      }
      
  out_loops:

      for(a=(cx+n_p_class-r)%n_p_class, aa=0;aa<t;a=(a+1)%n_p_class, aa++) 
      for(b=(cy+n_p_class-r)%n_p_class, bb=0;bb<t;b=(b+1)%n_p_class, bb++) {
         pointer = class_polar[tip][a][b];
         while(pointer != NULL) {
            comparisons ++ ;
            counter ++ ;
            s1 = pointer->sum;
            s2 = pointer->sum2;
            t1 = 0.0;
            i = pointer->ptr_x >> 1;
            j = pointer->ptr_y >> 1;
   
            switch(isom) {
             case IDENTITY   :  
                     for(x=0;x< size;x++)
                     for(y=0;y< size;y++)
                         t1 += contract[x+i][y+j] * range[x][y];
                     break;
             case R_ROTATE90 :  
                     for(x=0;x< size;x++)
                     for(y=0;y< size;y++)
                         t1 += contract[x+i][y+j] * range[y][size-x-1];
                     break;
             case L_ROTATE90 : 
                     for(x=0;x< size;x++)
                     for(y=0;y< size;y++)
                         t1 += contract[x+i][y+j] * range[size-y-1][x];
                     break;
             case ROTATE180  : 
                     for(x=0;x< size;x++)
                     for(y=0;y< size;y++)
                         t1 += contract[x+i][y+j] * range[size-x-1][size-y-1];
                     break;
             case R_VERTICAL : 
                     for(x=0;x< size;x++)
                     for(y=0;y< size;y++)
                         t1 += contract[x+i][y+j] * range[x][size-y-1];
                     break;
             case R_HORIZONTAL:
                     for(x=0;x< size;x++)
                     for(y=0;y< size;y++)
                         t1 += contract[x+i][y+j] * range[size-x-1][y];
                     break;
             case F_DIAGONAL:  
                     for(x=0;x< size;x++)
                     for(y=0;y< size;y++)
                         t1 += contract[x+i][y+j] * range[y][x];
                     break;
             case S_DIAGONAL:  
                     for(x=0;x< size;x++)
                     for(y=0;y< size;y++)
                         t1 += contract[x+i][y+j] * range[size-y-1][size-x-1];
                     break;
            }


            /* Compute the scalig factor */
            det = s2*s0 - s1*s1;
            if(det == 0.0)
               alfa = 0.0;
            else
              alfa = (s0*t1 - s1*t0) / det;
            if(alfa < 0.0 ) alfa = 0.0;

            /* Quantize the scalig factor */
            qalfa = 0.5 + (alfa )/( MAX_ALFA)* (1 << N_BITALFA);
            if(qalfa < 0)  qalfa = 0;
            if(qalfa >= (1 << N_BITALFA))  qalfa = (1 << N_BITALFA) - 1;

            /* Compute the scalig factor back from the quantized value*/
            alfa = (double) qalfa / (double)(1 << N_BITALFA) * ( MAX_ALFA) ;

            /* Compute the offset */
            beta= (t0 - alfa*s1) / s0;
            if (alfa > 0.0)  beta += alfa * 255;

            /* Quantize the offset */
            qbeta = 0.5 + beta/ ((1.0+fabs(alfa))*255)*((1 << N_BITBETA)-1);
            if (qbeta< 0) qbeta = 0;
            if (qbeta >= 1 << N_BITBETA) qbeta = (1 << N_BITBETA)-1;

            /* Compute the offset back from the quantized value*/
            beta = (double) qbeta/(double)((1 << N_BITBETA)-1)*
                                              ((1.0+fabs(alfa)) * 255);
            if (alfa > 0.0) beta  -= alfa * 255;

            /* Compute the rms distance */
            sum = t2 - 2*alfa*t1 -2*beta*t0 + alfa*alfa*s2 +
                              2*alfa*beta*s1 + s0*beta*beta;
            dist = sqrt(sum / s0);

            if(dist < min ) {
               min = dist;
               *xd = pointer->ptr_x;
               *yd = pointer->ptr_y;
               *is = isom;
               *qalf = qalfa;
               *qbet = qbeta;
            }
            pointer = pointer -> next;
         }
      }
  }
  return (min) ;
}



double SaupeCoding(int atx,int aty,int size,int *xd,int *yd,int *is,
                                                       int* qalf,int *qbet)
{
  int x,y,qalfa ,qbeta,i,j;
  int tip,ii,found,counter,isom;
  double dist,alfa,beta;
  double min = 1000000000.0;
  double sum,s0,mean,det;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  static float r_vector[4096]; 
  static int nlist[MAX_NEIGHBOURS]; 
  register double pixel;


  tip = (int) rint(log((double) size)/log (2.0));

  if(tip == 0) {   /* size = 1 */
       *qbet = (int)image[atx][aty];
       *qalf = zeroalfa;
       *xd = 0;
       *yd = 0;
       *is = IDENTITY;
       return 0.0;
  }

  s0 = size * size;
  for(x=0;x < size;x++)
  for(y=0;y < size;y++) {
     pixel = (double)image[atx+x][aty+y];
     t0 +=  pixel;
     t2 +=  pixel * pixel;
     range[x][y] = pixel;
  }

  mean = t0 / s0;

  for(isom=IDENTITY;isom <= S_DIAGONAL;isom++) {
     switch(isom) {     
         case L_ROTATE90  : flips(size,range,flip_range,R_ROTATE90);
                            break;
         case R_ROTATE90  : flips(size,range,flip_range,L_ROTATE90);
                            break;
         default          : flips(size,range,flip_range,isom);
                            break;

     }

     ComputeSaupeVectors(flip_range,size,tip,r_vector);

     found = kdtree_search(r_vector, f_vectors[tip], feat_vect_dim[tip], 
                                              kd_tree[tip], eps ,matches,nlist); 

     for(ii=0;ii<found ;ii++) {
        comparisons ++ ;
        counter ++ ;
        s1 = codebook[tip][nlist[ii]].sum;
        s2 = codebook[tip][nlist[ii]].sum2;
        t1 = 0.0;
        i = codebook[tip][nlist[ii]].ptr_x >> 1;
        j = codebook[tip][nlist[ii]].ptr_y >> 1;
   
        switch(isom) {
         case IDENTITY   :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j] * range[x][y];
                 break;
         case R_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j] * range[y][size-x-1];
                 break;
         case L_ROTATE90 :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j] * range[size-y-1][x];
                 break;
         case ROTATE180  :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j] * range[size-x-1][size-y-1];
                 break;
         case R_VERTICAL :  
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j] * range[x][size-y-1];
                 break;
         case R_HORIZONTAL: 
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j] * range[size-x-1][y];
                 break;
         case F_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j] * range[y][x];
                 break;
         case S_DIAGONAL:   
                 for(x=0;x< size;x++)
                 for(y=0;y< size;y++)
                     t1 += contract[x+i][y+j] * range[size-y-1][size-x-1];
                 break;
        }

        /* Compute the scalig factor */
        det = s2*s0 - s1*s1;
        if(det == 0.0)
           alfa = 0.0;
        else
          alfa = (s0*t1 - s1*t0) / det;
        if(alfa < 0.0 ) alfa = 0.0;

        /* Quantize the scalig factor */
        qalfa = 0.5 + (alfa )/( MAX_ALFA)* (1 << N_BITALFA);
        if(qalfa < 0)  qalfa = 0;
        if(qalfa >= (1 << N_BITALFA))  qalfa = (1 << N_BITALFA) - 1;

        /* Compute the scalig factor back from the quantized value*/
        alfa = (double) qalfa / (double)(1 << N_BITALFA) * ( MAX_ALFA) ;

        /* Compute the offset */
        beta= (t0 - alfa*s1) / s0;
        if (alfa > 0.0)  beta += alfa * 255;

        /* Quantize the offset */
        qbeta = 0.5 + beta/ ((1.0+fabs(alfa))*255)*((1 << N_BITBETA)-1);
        if (qbeta< 0) qbeta = 0;
        if (qbeta >= 1 << N_BITBETA) qbeta = (1 << N_BITBETA)-1;

        /* Compute the offset back from the quantized value*/
        beta = (double) qbeta/(double)((1 << N_BITBETA)-1)*
                                          ((1.0+fabs(alfa)) * 255);
        if (alfa > 0.0) beta  -= alfa * 255;

        /* Compute the rms distance */
        sum = t2 - 2*alfa*t1 -2*beta*t0 + alfa*alfa*s2 +
                              2*alfa*beta*s1 + s0*beta*beta;
        dist = sqrt(sum / s0);

        if(dist < min ) {
            min = dist;
            *xd = codebook[tip][nlist[ii]].ptr_x;
            *yd = codebook[tip][nlist[ii]].ptr_y; 
            *is = isom;
            *qalf = qalfa;
            *qbet = qbeta;
        }
     }
  }
  return (min) ;
}

int best_rmean(int ics,int yps,int size,double alfa)
{
  double rmean;
  double sum = 0.0;
  int qrmean;
  int x,y;

  for(x=0;x < size;x++)
  for(y=0;y < size;y++) 
     sum += (double)image[ics+x][yps+y];

  rmean = sum / (double)(size * size);
   qrmean = quan(rmean);


  return qrmean;
}

int best_beta(int ics,int yps,int size,double alfa)
{
  double beta;
  double sum = 0.0;
  int qbeta;
  int x,y;

  for(x=0;x < size;x++)
  for(y=0;y < size;y++) 
     sum += (double)image[ics+x][yps+y];

  beta = sum / (double)(size * size);
  if (alfa > 0.0)  beta += alfa * 255;

  qbeta = 0.5 + beta/ ((1.0+fabs(alfa))*255)*((1 << N_BITBETA)-1);
  if (qbeta< 0) qbeta = 0;
  if (qbeta >= 1 << N_BITBETA) qbeta = (1 << N_BITBETA)-1;

  return qbeta;
}



void testing_quadtree(int atx,int aty,int size,double tol_entr,
                                       double tol_rms,double tol_var)
{
  double best_rms;
  int domx,domy,isom,qalfa,qbeta,qrmean;
  static int coded = 0;
  static int oldper = 0;
  static int newper = 0;
  double compress,bpp;
  int bytes;
  int k;

  if(atx >= image_height  || aty >= image_width )
      return;

  if (size > max_size || atx+size > image_height || aty+size > image_width ) {

     for(k=aty;k<aty+size;k++)
       qtt[atx+size/2][k] = 0;

     for(k=atx;k<atx+size;k++)
       qtt[k][aty+size/2] = 0;

     testing_quadtree(atx,aty,size/2, tol_entr,tol_rms,tol_var);
     testing_quadtree(atx+size/2,aty,size/2,tol_entr,tol_rms,tol_var);
     testing_quadtree(atx,aty+size/2,size/2,tol_entr,tol_rms,tol_var);
     testing_quadtree(atx+size/2,aty+size/2,size/2,tol_entr,tol_rms,tol_var);
     return;
  }

  if ((entropy(size,size, atx, aty) > tol_entr ||
                   variance(size,size, atx, aty) > tol_var) && size > min_size) {

     tol_entr = tol_entr + (log((double)adapt)/log(2.0)) / 
               (log((double)max_size)/log(2.0) - log((double)size)/log(2.0) + 1) ;

     pack(1,(long)1,fp);
     for(k=aty;k<aty+size;k++)
        qtt[atx+size/2][k] = 0;

     for(k=atx;k<atx+size;k++)
         qtt[k][aty+size/2] = 0;   

     testing_quadtree(atx,aty,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
     testing_quadtree(atx+size/2,aty,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
     testing_quadtree(atx,aty+size/2,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
     testing_quadtree(atx+size/2,aty+size/2,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);

  } else {

     best_rms = testing_Coding(atx,aty,size,&domx,&domy,&isom,&qalfa,&qrmean);
     if (best_rms > tol_rms  && size > min_size) {

        tol_entr =  tol_entr + (log(adapt)/log(2.0)) /
                 (log((double)max_size)/log(2.0) - log((double)size)/log(2.0) + 1) ;

        pack(1,(long)1,fp);
        for(k=aty;k<aty+size;k++)
           qtt[atx+size/2][k] = 0;

        for(k=atx;k<atx+size;k++)
           qtt[k][aty+size/2] = 0;

        testing_quadtree(atx,aty,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
        testing_quadtree(atx+size/2,aty,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
        testing_quadtree(atx,aty+size/2,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
        testing_quadtree(atx+size/2,aty+size/2,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);

     } else {

        if(size > min_size)
            pack(1,(long)0,fp);

        if(abs(qalfa - zeroalfa) <= zero_threshold) {
            qrmean = best_beta(atx,aty,size,0.0);
            qalfa = zeroalfa;
        }
        pack(N_BITALFA, (long)qalfa, fp);
        pack(N_BITRMEAN, (long)qrmean
          , fp);
        if(isColor){
          double um = u_mean(size,size,atx,aty);
          double vm = v_mean(size,size,atx,aty);
          pack(8,(long)um,fp);
          pack(8,(long)vm,fp);
        }
        // printf("qalfa = %d\n",qalfa);
        zero_alfa_transform ++;
        if(qalfa != zeroalfa) {
           zero_alfa_transform --;
           pack(3, (long)isom, fp);
           pack(bits_per_coordinate_h,(long)(domx / SHIFT),fp);
           pack(bits_per_coordinate_w,(long)(domy / SHIFT),fp);
        }

        transforms ++;
        coded += size * size ;
        newper = (int) (((double) coded / (image_width * image_height)) * 100);
        if(newper > oldper) {
           bytes = pack(-2,(long)0,fp);
           compress =(double) coded / (double)bytes;
           bpp      = 8.0 /compress;
           printf(" %d %% Coded, Compression %f:1, bpp %f\r", 
                                      newper,compress,bpp);
           oldper = newper;
           fflush(stdout);
        }
     }
  }
}

void LumInv_quadtree(int atx,int aty,int size,double tol_entr,
                                       double tol_rms,double tol_var)
{
  double best_rms;
  int domx,domy,isom,qalfa,qrmean;
  static int coded = 0;
  static int oldper = 0;
  static int newper = 0;
  double compress,bpp;
  int bytes;
  int k;

  if(atx >= image_height  || aty >= image_width )
      return;

  if (size > max_size || atx+size > image_height || aty+size > image_width ) {

     for(k=aty;k<aty+size;k++)
       qtt[atx+size/2][k] = 0;

     for(k=atx;k<atx+size;k++)
       qtt[k][aty+size/2] = 0;

     LumInv_quadtree(atx,aty,size/2, tol_entr,tol_rms,tol_var);
     LumInv_quadtree(atx+size/2,aty,size/2,tol_entr,tol_rms,tol_var);
     LumInv_quadtree(atx,aty+size/2,size/2,tol_entr,tol_rms,tol_var);
     LumInv_quadtree(atx+size/2,aty+size/2,size/2,tol_entr,tol_rms,tol_var);
     return;
  }

  if ((entropy(size,size, atx, aty) > tol_entr ||
                   variance(size,size, atx, aty) > tol_var) && size > min_size) {

     tol_entr = tol_entr + (log((double)adapt)/log(2.0)) / 
               (log((double)max_size)/log(2.0) - log((double)size)/log(2.0) + 1) ;

     pack(1,(long)1,fp);
     for(k=aty;k<aty+size;k++)
        qtt[atx+size/2][k] = 0;

     for(k=atx;k<atx+size;k++)
         qtt[k][aty+size/2] = 0;   

     LumInv_quadtree(atx,aty,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
     LumInv_quadtree(atx+size/2,aty,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
     LumInv_quadtree(atx,aty+size/2,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
     LumInv_quadtree(atx+size/2,aty+size/2,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);

  } else {

     best_rms = LumInv_Coding(atx,aty,size,&domx,&domy,&isom,&qalfa,&qrmean);
     if (best_rms > tol_rms  && size > min_size) {

        tol_entr =  tol_entr + (log(adapt)/log(2.0)) /
                 (log((double)max_size)/log(2.0) - log((double)size)/log(2.0) + 1) ;

        pack(1,(long)1,fp);
        for(k=aty;k<aty+size;k++)
           qtt[atx+size/2][k] = 0;

        for(k=atx;k<atx+size;k++)
           qtt[k][aty+size/2] = 0;

        LumInv_quadtree(atx,aty,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
        LumInv_quadtree(atx+size/2,aty,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
        LumInv_quadtree(atx,aty+size/2,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
        LumInv_quadtree(atx+size/2,aty+size/2,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);

     } else {

        if(size > min_size)
            pack(1,(long)0,fp);

        if(abs(qalfa - zeroalfa) <= zero_threshold) {
            
            qalfa = zeroalfa;
            // qrmean = best_rmean(atx,aty,size,0.0);
        }
        pack(N_BITALFA, (long)qalfa, fp);
        pack(N_BITRMEAN, (long)qrmean, fp);
        if(isColor){
          double um = u_mean(size,size,atx,aty);
          double vm = v_mean(size,size,atx,aty);
          pack(8,(long)um,fp);
          pack(8,(long)vm,fp);
        }
        zero_alfa_transform ++;
        // printf("qalfa = %d\n",qalfa);
        if(qalfa != zeroalfa) {
           zero_alfa_transform --;
           pack(3, (long)isom, fp);
           pack(bits_per_coordinate_h,(long)(domx / SHIFT),fp);
           pack(bits_per_coordinate_w,(long)(domy / SHIFT),fp);
        }

        transforms ++;
        coded += size * size ;
        newper = (int) (((double) coded / (image_width * image_height)) * 100);
        if(newper > oldper) {
           bytes = pack(-2,(long)0,fp);
           compress =(double) coded / (double)bytes;
           bpp      = 8.0 /compress;
           printf(" %d %% Coded, Compression %f:1, bpp %f\r", 
                                      newper,compress,bpp);
           oldper = newper;
           fflush(stdout);
        }
     }
  }
}

void Nonlinear_quadtree(int atx,int aty,int size,double tol_entr,
                                       double tol_rms,double tol_var)
{
  double best_rms;
  int domx,domy,isom,qalfa1,qalfa2,qbeta;
  static int coded = 0;
  static int oldper = 0;
  static int newper = 0;
  double compress,bpp;
  int bytes;
  int k;

  if(atx >= image_height  || aty >= image_width )
      return;

  if (size > max_size || atx+size > image_height || aty+size > image_width ) {

     for(k=aty;k<aty+size;k++)
       qtt[atx+size/2][k] = 0;

     for(k=atx;k<atx+size;k++)
       qtt[k][aty+size/2] = 0;

     Nonlinear_quadtree(atx,aty,size/2, tol_entr,tol_rms,tol_var);
     Nonlinear_quadtree(atx+size/2,aty,size/2,tol_entr,tol_rms,tol_var);
     Nonlinear_quadtree(atx,aty+size/2,size/2,tol_entr,tol_rms,tol_var);
     Nonlinear_quadtree(atx+size/2,aty+size/2,size/2,tol_entr,tol_rms,tol_var);
     return;
  }

  if ((entropy(size,size, atx, aty) > tol_entr ||
                   variance(size,size, atx, aty) > tol_var) && size > min_size) {

     tol_entr = tol_entr + (log((double)adapt)/log(2.0)) / 
               (log((double)max_size)/log(2.0) - log((double)size)/log(2.0) + 1) ;

     pack(1,(long)1,fp);
     for(k=aty;k<aty+size;k++)
        qtt[atx+size/2][k] = 0;

     for(k=atx;k<atx+size;k++)
         qtt[k][aty+size/2] = 0;   

     Nonlinear_quadtree(atx,aty,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
     Nonlinear_quadtree(atx+size/2,aty,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
     Nonlinear_quadtree(atx,aty+size/2,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
     Nonlinear_quadtree(atx+size/2,aty+size/2,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);

  } else {

     best_rms = Nonlinear_Coding(atx,aty,size,&domx,&domy,&isom,&qalfa1,&qalfa2,&qbeta);
     if (best_rms > tol_rms  && size > min_size) {

        tol_entr =  tol_entr + (log(adapt)/log(2.0)) /
                 (log((double)max_size)/log(2.0) - log((double)size)/log(2.0) + 1) ;

        pack(1,(long)1,fp);
        for(k=aty;k<aty+size;k++)
           qtt[atx+size/2][k] = 0;

        for(k=atx;k<atx+size;k++)
           qtt[k][aty+size/2] = 0;

        Nonlinear_quadtree(atx,aty,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
        Nonlinear_quadtree(atx+size/2,aty,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
        Nonlinear_quadtree(atx,aty+size/2,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
        Nonlinear_quadtree(atx+size/2,aty+size/2,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);

     } else {

        if(size > min_size)
            pack(1,(long)0,fp);

        // if(abs(qalfa1 - zeroalfa) <= zero_threshold) {
        //     qbeta = best_beta(atx,aty,size,0.0);
        //     qalfa1 = zeroalfa;
        // }

        pack(N_BITALFA1, (long)qalfa1, fp);
        pack(N_BITALFA2, (long)qalfa2, fp);
        pack(N_BITBETA2, (long)qbeta, fp);
        if(isColor){
          double um = u_mean(size,size,atx,aty);
          double vm = v_mean(size,size,atx,aty);
          pack(8,(long)um,fp);
          pack(8,(long)vm,fp);
        }
        zero_alfa_transform ++;
        if(qalfa1 != zeroalfa || qalfa2 != zeroalfa) {
           zero_alfa_transform --;
           pack(3, (long)isom, fp);
           pack(bits_per_coordinate_h,(long)(domx / SHIFT),fp);
           pack(bits_per_coordinate_w,(long)(domy / SHIFT),fp);
        }

        transforms ++;
        coded += size * size ;
        newper = (int) (((double) coded / (image_width * image_height)) * 100);
        if(newper > oldper) {
           bytes = pack(-2,(long)0,fp);
           compress =(double) coded / (double)bytes;
           bpp      = 8.0 /compress;
           printf(" %d %% Coded, Compression %f:1, bpp %f\r", 
                                      newper,compress,bpp);
           oldper = newper;
           fflush(stdout);
        }
     }
  }
}


int  bitlength(unsigned long val)
{
    int bits = 1;
    
    if (val > 0xffff) bits += 16, val >>= 16;
    if (val > 0xff)   bits += 8,  val >>= 8;
    if (val > 0xf)    bits += 4,  val >>= 4;
    if (val > 0x3)    bits += 2,  val >>= 2;
    if (val > 0x1)    bits += 1;
    return bits;
}

/* paritioning based on Range size */
void quadtree_2(int atx,int aty,int size,double tol_entr,
                                       double tol_rms,double tol_var)
{
  double best_rms;
  int domx,domy,isom,qalfa,qbeta;
  static int coded = 0;
  static int oldper = 0;
  static int newper = 0;
  double compress,bpp;
  int bytes;
  int k;

 /* the is modification added to quadtree*/
  int s_log, s_size;
  s_log = bitlength(size) - 1;
  s_size = 1 << s_log;
  // printf("s_size = %d\n",s_size);

  if(atx >= image_height  || aty >= image_width )
      return;

  if (size > max_size || atx+size > image_height || aty+size > image_width ) {

     for(k=aty;k<aty+size;k++)
       qtt[atx+size/2][k] = 0;

     for(k=atx;k<atx+size;k++)
       qtt[k][aty+size/2] = 0;

     quadtree_2(atx,aty,size/2, tol_entr,tol_rms,tol_var);
     quadtree_2(atx+size/2,aty,size/2,tol_entr,tol_rms,tol_var);
     quadtree_2(atx,aty+size/2,size/2,tol_entr,tol_rms,tol_var);
     quadtree_2(atx+size/2,aty+size/2,size/2,tol_entr,tol_rms,tol_var);
     return;
  }

  if ((entropy(size,size, atx, aty) > tol_entr ||
                   variance(size,size, atx, aty) > tol_var) && size > min_size) {

     tol_entr = tol_entr + (log((double)adapt)/log(2.0)) / 
               (log((double)max_size)/log(2.0) - log((double)size)/log(2.0) + 1) ;

     pack(1,(long)1,fp);
     for(k=aty;k<aty+size;k++)
        qtt[atx+size/2][k] = 0;

     for(k=atx;k<atx+size;k++)
         qtt[k][aty+size/2] = 0;   

     quadtree_2(atx,aty,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
     quadtree_2(atx+size/2,aty,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
     quadtree_2(atx,aty+size/2,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
     quadtree_2(atx+size/2,aty+size/2,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
  }else if (s_log > MAX_BITS){

     pack(1,(long)1,fp);
     for(k=aty;k<aty+size;k++)
        qtt[atx+size/2][k] = 0;

     for(k=atx;k<atx+size;k++)
         qtt[k][aty+size/2] = 0;   

     quadtree_2(atx,aty,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
     quadtree_2(atx+size/2,aty,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
     quadtree_2(atx,aty+size/2,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
     quadtree_2(atx+size/2,aty+size/2,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);

  } else {

     best_rms = Coding(atx,aty,size,&domx,&domy,&isom,&qalfa,&qbeta);
     if (best_rms > tol_rms  && size > min_size) {

        tol_entr =  tol_entr + (log(adapt)/log(2.0)) /
                 (log((double)max_size)/log(2.0) - log((double)size)/log(2.0) + 1) ;

        pack(1,(long)1,fp);
        for(k=aty;k<aty+size;k++)
           qtt[atx+size/2][k] = 0;

        for(k=atx;k<atx+size;k++)
           qtt[k][aty+size/2] = 0;

        quadtree_2(atx,aty,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
        quadtree_2(atx+size/2,aty,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
        quadtree_2(atx,aty+size/2,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
        quadtree_2(atx+size/2,aty+size/2,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);

     } else {

        if(size > min_size)
            pack(1,(long)0,fp);

        if(abs(qalfa - zeroalfa) <= zero_threshold) {
            qbeta = best_beta(atx,aty,size,0.0);
            qalfa = zeroalfa;
        }
        pack(N_BITALFA, (long)qalfa, fp);
        pack(N_BITBETA, (long)qbeta, fp);
        if(isColor){
          double um = u_mean(size,size,atx,aty);
          double vm = v_mean(size,size,atx,aty);
          pack(8,(long)um,fp);
          pack(8,(long)vm,fp);
        }
        // printf("qalfa = %d\n",qalfa);
        zero_alfa_transform ++;
        if(qalfa != zeroalfa) {
           zero_alfa_transform --;
           pack(3, (long)isom, fp);
           pack(bits_per_coordinate_h,(long)(domx / SHIFT),fp);
           pack(bits_per_coordinate_w,(long)(domy / SHIFT),fp);
        }

        transforms ++;
        coded += size * size ;
        newper = (int) (((double) coded / (image_width * image_height)) * 100);
        if(newper > oldper) {
           bytes = pack(-2,(long)0,fp);
           compress =(double) coded / (double)bytes;
           bpp      = 8.0 /compress;
           printf(" %d %% Coded, Compression %f:1, bpp %f\r", 
                                      newper,compress,bpp);
           COMPRESS = compress;
           BPP = bpp;
           oldper = newper;
           fflush(stdout);
        }
     }
  }
}

void HV_patition(int atx,int aty,int size)
{
  int s_log = (int)log2(size);
  int s_size = 1 << s_log;
  double best_rms;
  int domx,domy,isom,qalfa,qbeta;
  static int coded = 0;
  static int oldper = 0;
  static int newper = 0;
  double compress,bpp;
  int bytes;
  int k;

  
  if(atx >= image_height  || aty >= image_width )
      return;

  if (size > max_size || atx+size > image_height || aty+size > image_width ) {
     
     HV_patition(atx,aty,size/2);
     HV_patition(atx+size/2,aty,size/2);
     HV_patition(atx,aty+size/2,size/2);
     HV_patition(atx+size/2,aty+size/2,size/2);
     return;

  } else {

     best_rms = Coding(atx,aty,size,&domx,&domy,&isom,&qalfa,&qbeta);
     if (size > min_size) {
      pack(1,(long)0,fp); // encode affine map

      if(abs(qalfa - zeroalfa) <= zero_threshold) {
               qbeta = best_beta(atx,aty,size,0.0);
               qalfa = zeroalfa;
            }
            pack(N_BITALFA, (long)qalfa, fp);
            pack(N_BITBETA, (long)qbeta, fp);
            if(isColor){
               double um = u_mean(size,size,atx,aty);
               double vm = v_mean(size,size,atx,aty);
               pack(8,(long)um,fp);
               pack(8,(long)vm,fp);
            }
        // printf("qalfa = %d\n",qalfa);
            zero_alfa_transform ++;
            if(qalfa != zeroalfa) {
               zero_alfa_transform --;
               pack(3, (long)isom, fp);
               pack(bits_per_coordinate_h,(long)(domx / SHIFT),fp);
               pack(bits_per_coordinate_w,(long)(domy / SHIFT),fp);
            }

           transforms ++;
           coded += size * size ;
           newper = (int) (((double) coded / (image_width * image_height)) * 100);
           if(newper > oldper) {
              bytes = pack(-2,(long)0,fp);
              compress =(double) coded / (double)bytes;
              bpp      = 8.0 /compress;
              printf(" %d %% Coded, Compression %f:1, bpp %f\r", 
			                                newper,compress,bpp);
              COMPRESS = compress;
              BPP = bpp;
              oldper = newper;
              fflush(stdout);
           }
    } else {

           if(size > min_size)
               pack(1,(long)1,fp);
           
           int x1,y1,x2,y2,Rp,Cp;
	   double h,v,hmax,vmax;
	   hmax = vmax = 0; 
           
    }
  }
}




void quadtree(int atx,int aty,int size,double tol_entr,
                                       double tol_rms,double tol_var)
{
  double best_rms;
  int domx,domy,isom,qalfa,qbeta;
  static int coded = 0;
  static int oldper = 0;
  static int newper = 0;
  double compress,bpp;
  int bytes;
  int k;

  if(atx >= image_height  || aty >= image_width )
      return;

  if (size > max_size || atx+size > image_height || aty+size > image_width ) {

     for(k=aty;k<aty+size;k++)
       qtt[atx+size/2][k] = 0;

     for(k=atx;k<atx+size;k++)
       qtt[k][aty+size/2] = 0;

     quadtree(atx,aty,size/2, tol_entr,tol_rms,tol_var);
     quadtree(atx+size/2,aty,size/2,tol_entr,tol_rms,tol_var);
     quadtree(atx,aty+size/2,size/2,tol_entr,tol_rms,tol_var);
     quadtree(atx+size/2,aty+size/2,size/2,tol_entr,tol_rms,tol_var);
     return;
  }

  if ((entropy(size,size, atx, aty) > tol_entr ||
                   variance(size,size, atx, aty) > tol_var) && size > min_size) {

     tol_entr = tol_entr + (log((double)adapt)/log(2.0)) / 
               (log((double)max_size)/log(2.0) - log((double)size)/log(2.0) + 1) ;

     pack(1,(long)1,fp);
     for(k=aty;k<aty+size;k++)
        qtt[atx+size/2][k] = 0;

     for(k=atx;k<atx+size;k++)
         qtt[k][aty+size/2] = 0;   

     quadtree(atx,aty,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
     quadtree(atx+size/2,aty,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
     quadtree(atx,aty+size/2,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
     quadtree(atx+size/2,aty+size/2,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);

  } else {

     best_rms = Coding(atx,aty,size,&domx,&domy,&isom,&qalfa,&qbeta);
     if (best_rms > tol_rms  && size > min_size) {

        tol_entr =  tol_entr + (log(adapt)/log(2.0)) /
                 (log((double)max_size)/log(2.0) - log((double)size)/log(2.0) + 1) ;

        pack(1,(long)1,fp);
        for(k=aty;k<aty+size;k++)
           qtt[atx+size/2][k] = 0;

        for(k=atx;k<atx+size;k++)
           qtt[k][aty+size/2] = 0;

        quadtree(atx,aty,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
        quadtree(atx+size/2,aty,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
        quadtree(atx,aty+size/2,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
        quadtree(atx+size/2,aty+size/2,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);

     } else {

        if(size > min_size)
            pack(1,(long)0,fp);

        if(abs(qalfa - zeroalfa) <= zero_threshold) {
            qbeta = best_beta(atx,aty,size,0.0);
            qalfa = zeroalfa;
        }
        pack(N_BITALFA, (long)qalfa, fp);
        pack(N_BITBETA, (long)qbeta, fp);
        if(isColor){
          double um = u_mean(size,size,atx,aty);
          double vm = v_mean(size,size,atx,aty);
          pack(8,(long)um,fp);
          pack(8,(long)vm,fp);
        }
        // printf("qalfa = %d\n",qalfa);
        zero_alfa_transform ++;
        if(qalfa != zeroalfa) {
           zero_alfa_transform --;
           pack(3, (long)isom, fp);
           pack(bits_per_coordinate_h,(long)(domx / SHIFT),fp);
           pack(bits_per_coordinate_w,(long)(domy / SHIFT),fp);
        }

        transforms ++;
        coded += size * size;
        newper = (int) (((double) coded / (image_width * image_height)) * 100);
        if(newper > oldper) {
           bytes = pack(-2,(long)0,fp);
           compress =(double) coded / (double)bytes;
           bpp      = 8.0 /compress;
           printf(" %d %% Coded, Compression %f:1, bpp %f\r", 
			                                newper,compress,bpp);
           COMPRESS = compress;
           BPP = bpp;
           oldper = newper;
           fflush(stdout);
        }
     }
  }
}

void adaptiveFlips(int x_size,int y_size,double **block,double **flip_block,int iso)
{
  register int i,j;

  switch(iso) {
     case IDENTITY   :  
             for(i=0;i< y_size;i++)
             for(j=0;j< x_size;j++)
                 flip_block[i][j] = block[i][j];
             break;
     case L_ROTATE90 :  
             for(i=0;i< y_size;i++)
             for(j=0;j< x_size;j++)
                 flip_block[i][j] = block[j][x_size-i-1];
             break;
     case R_ROTATE90 :  
             for(i=0;i< y_size;i++)
             for(j=0;j< x_size;j++)
                 flip_block[i][j] = block[y_size-j-1][i];
             break;
     case ROTATE180  :  
             for(i=0;i< y_size;i++)
             for(j=0;j< x_size;j++)
                 flip_block[i][j] = block[y_size-i-1][x_size-j-1];
             break;
     case R_VERTICAL :  
             for(i=0;i< y_size;i++)
             for(j=0;j< x_size;j++)
                 flip_block[i][j] = block[i][x_size-j-1];
             break;
     case R_HORIZONTAL: 
             for(i=0;i< y_size;i++)
             for(j=0;j< x_size;j++)
                 flip_block[i][j] = block[y_size-i-1][j];
             break;
     case F_DIAGONAL:   
             for(i=0;i< y_size;i++)
             for(j=0;j< x_size;j++)
                 flip_block[i][j] = block[j][i];
             break;
     case S_DIAGONAL:   
             for(i=0;i< y_size;i++)
             for(j=0;j< x_size;j++)
                 flip_block[i][j] = block[y_size-j-1][x_size-i-1];
             break;
   }
}

void flips(int size,double **block,double **flip_block,int iso)
{
  register int i,j;

  switch(iso) {
     case IDENTITY   :  
             for(i=0;i< size;i++)
             for(j=0;j< size;j++)
                 flip_block[i][j] = block[i][j];
             break;
     case L_ROTATE90 :  
             for(i=0;i< size;i++)
             for(j=0;j< size;j++)
                 flip_block[i][j] = block[j][size-i-1];
             break;
     case R_ROTATE90 :  
             for(i=0;i< size;i++)
             for(j=0;j< size;j++)
                 flip_block[i][j] = block[size-j-1][i];
             break;
     case ROTATE180  :  
             for(i=0;i< size;i++)
             for(j=0;j< size;j++)
                 flip_block[i][j] = block[size-i-1][size-j-1];
             break;
     case R_VERTICAL :  
             for(i=0;i< size;i++)
             for(j=0;j< size;j++)
                 flip_block[i][j] = block[i][size-j-1];
             break;
     case R_HORIZONTAL: 
             for(i=0;i< size;i++)
             for(j=0;j< size;j++)
                 flip_block[i][j] = block[size-i-1][j];
             break;
     case F_DIAGONAL:   
             for(i=0;i< size;i++)
             for(j=0;j< size;j++)
                 flip_block[i][j] = block[j][i];
             break;
     case S_DIAGONAL:   
             for(i=0;i< size;i++)
             for(j=0;j< size;j++)
                 flip_block[i][j] = block[size-j-1][size-i-1];
             break;
   }
}





void traverseImage(int atx, int aty, int x_size, int y_size)
{
  
  int s_size;
  int s_log;

  s_log = (int) log2(min_2(x_size,y_size));
  s_size = 1 << s_log;
  // printf("s_size = %d\n",s_size);

  // printf("s_log = %d\n",s_log);
  if(s_log > MAX_ADAP_R_BITS){
    traverseImage(atx,          aty,          s_size/2, s_size/2);
    traverseImage(atx+s_size/2, aty,          s_size/2, s_size/2);
    traverseImage(atx,          aty+s_size/2, s_size/2, s_size/2);
    traverseImage(atx+s_size/2, aty+s_size/2, s_size/2, s_size/2);
  }
  else{
    compressRange(atx,aty, x_size, y_size);
  }

  if(x_size > s_size)
    traverseImage(atx+s_size,aty,x_size - s_size, y_size);
  
  if(y_size > s_size)
    traverseImage(atx, aty + s_size,s_size, y_size - s_size);
}


void compressRange(int atx, int aty, int x_size, int y_size)
{
  if(atx >= image_height  || aty >= image_width )
      return;


  if(x_size == 0 || y_size == 0)
    return; 


  // printf("isnide compress range:\n");
  // printf("x_size = %d \t y_size = %d\n",x_size,y_size);
  int domx,domy,isom,qalfa,qbeta;
  static int coded = 0;
  static int oldper = 0;
  static int newper = 0;
  static int else_count = 0;
  double compress,bpp;
  int bytes;
  int k;
  // printf("X size = %d \t Y size = %d \n",x_size,y_size);
  

  double best_rms;
  best_rms = adaptiveFisherCoding(atx,aty,x_size,y_size,&domx,&domy,&isom,&qalfa,&qbeta);
  // printf("best rms = %f\n",best_rms);

  // exit(0);
  if(max_2(x_size,y_size) == 1 << (MIN_ADAP_R_BITS)  || best_rms < max_error2 * (double)(x_size * y_size)){
    // printf("packing bits \n");
    // printf("transforms = %d \n",transforms);
    // printf("x_sizse = %d \t y_size = %d\n",x_size,y_size);
    if(max_2(x_size,y_size) !=  1 << MIN_ADAP_R_BITS)
      pack(1,(long)1,fp);

    // if(abs(qalfa - zeroalfa) <= zero_threshold) {
    //     qbeta = best_beta(atx,aty,x_size,0.0);
    //     qalfa = zeroalfa;
    // }
    pack(N_BITALFA, (long)qalfa, fp);
    pack(N_BITBETA, (long)qbeta, fp);
    // if(isColor){
    //   double um = u_mean(x_size,y_size,atx,aty);
    //   double vm = v_mean(x_size,y_size,atx,aty);
    //   pack(8,(long)um,fp);
    //   pack(8,(long)vm,fp);
    // }
    // printf("qalfa = %d\n",qalfa);
    zero_alfa_transform ++;
    if(qalfa != zeroalfa) {
       zero_alfa_transform --;
       pack(3, (long)isom, fp);
       pack(bits_per_coordinate_h,(long)(domx / SHIFT),fp);
       pack(bits_per_coordinate_w,(long)(domy / SHIFT),fp);
    }

    transforms ++;
    // printf("transforms = %d\r",transforms);
    coded += x_size * y_size ;
    newper = (int) (((double) coded / (image_width * image_height)) * 100);
    if(newper > oldper) {
       bytes = pack(-2,(long)0,fp);
       compress =(double) coded / (double)bytes;
       bpp      = 8.0 /compress;
       printf(" %d %% Coded, Compression %f:1, bpp %f\r", 
                                  newper,compress,bpp);
       oldper = newper;
       fflush(stdout);
    }
  }
  else{
      // else_count++;
      // if(else_count == 1)
      // printf("x size = %d \t y size = %d\n",x_size,y_size);

    // if(x_size == 1 || y_size == 1){ 
    //   printf("x size = %d \t y size = %d\n",x_size,y_size);
    //    exit(0);
    // }  
    pack(1,(long)0,fp);
    int i,j,Rp, Cp1, Cp2, x1,x2,x3,x4,y1,y2 ,vmax1,vmax2,hmax,h,v1,v2;
    
    hmax=vmax1=vmax2=0;
 

    Cp1=Cp2=atx;
    Rp=aty;
    for(j=aty;j<aty+y_size-1;j++)
    {
          y1=y2=0;
          for(i=atx;i<atx+x_size;i++)
          {
                y1=y1+image[j][i];
                y2=y2+image[j+1][i];
          }
          h=(min_2(j-aty+1,y_size-j+aty-1)*abs(y1-y2));
          if(h>hmax)
          {
                hmax=h;         
                Rp=j;
          }
    }
    for(i=atx;i<atx+x_size-1;i++)
    {
          x1=x2=0;
          for(j=aty;j<aty+y_size;j++)
          {
                x1=x1+image[j][i];
                x2=x2+image[j][i+1];
          }
          v1=(min_2(i-atx+1,x_size-i+atx-1)*abs(x1-x2));
          if(v1>vmax1)
          {
                vmax1=v1;
                Cp1=i;
          }
          // x3=x4=0;
          // for(j=Rp;j<aty+y_size;j++)
          // {
          //       x3=x3+image[j][i];
          //       x3=x4+image[j][i+1];
          // }
          // v2=(min_2(i-atx+1,x_size-i+atx-1)*abs(x3-x4));
          // if(v2>vmax2)
          // {
          //       vmax2=v2;
          //       Cp2=i;
          // }
    }
   
      pack(3,(long)Rp-aty+1,fp);
      pack(3,(long)Cp1-atx+1,fp);
   //   pack(3,(long)Cp2-atx+1,fp);
     
      // printf("y_size1 =%d\n",Rp-aty+1,fp);
      // printf("x_size1=%d\n",Cp1-atx+1,fp);
      // printf("",Cp2-atx+1,fp);
     

      compressRange(atx,   aty,   Cp1-atx+1,        Rp-aty+1);
      compressRange(Cp1+1, aty,   atx+x_size-Cp1-1, Rp-aty+1);
      compressRange(atx,   Rp+1,  Cp1-atx+1,        aty+y_size-Rp-1);
      compressRange(Cp1+1, Rp+1,  atx+x_size-Cp1-1, aty+y_size-Rp-1);
    
  }
}




// void adaptive_quadtree(int atx,int aty,int x_size,int y_size)
// {
//   double best_rms;
//   int domx,domy,isom,qalfa,qbeta;
//   static int coded = 0;
//   static int oldper = 0;
//   static int newper = 0;
//   double compress,bpp;
//   int bytes;
//   int k;
//    int s_log = (int)log2(min_2(x_size,y_size));
//    int s_size = 1 << s_log;
//   // int s_size;
//   // int s_log;
//   // s_log = (int) log2(size);
//   // s_size = 1 << s_log;

//   if(atx >= image_height  || aty >= image_width || s_size == 0)
//       return;

//   if (s_size > max_size || atx+s_size > image_height || aty+s_size > image_width ) {

//      // for(k=aty;k<aty+size;k++)
//      //   qtt[atx+size/2][k] = 0;

//      // for(k=atx;k<atx+size;k++)
//      //   qtt[k][aty+size/2] = 0;

//      adaptive_quadtree(atx,aty,s_size/2,s_size/2);
//      adaptive_quadtree(atx+s_size/2,aty,s_size/2,s_size/2);
//      adaptive_quadtree(atx,aty+s_size/2,s_size/2,s_size/2);
//      adaptive_quadtree(atx+s_size/2,aty+s_size/2,s_size/2,s_size/2);
//      return;
//   //}

//   // if ((entropy(size,size, atx, aty) > tol_entr ||
//   //                  variance(size,size, atx, aty) > tol_var) && size > min_size) {

//   //    tol_entr = tol_entr + (log((double)adapt)/log(2.0)) / 
//   //              (log((double)max_size)/log(2.0) - log((double)size)/log(2.0) + 1) ;

//   //    pack(1,(long)1,fp);
//   //    for(k=aty;k<aty+size;k++)
//   //       qtt[atx+size/2][k] = 0;

//   //    for(k=atx;k<atx+size;k++)
//   //        qtt[k][aty+size/2] = 0;   

//   //    adaptive_quadtree(atx,aty,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
//   //    adaptive_quadtree(atx+size/2,aty,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
//   //    adaptive_quadtree(atx,aty+size/2,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);
//   //    adaptive_quadtree(atx+size/2,aty+size/2,size/2, tol_entr,tol_rms*adapt,tol_var*adapt);

//   } else {

//      best_rms = Coding(atx,aty,x_size,y_size,&domx,&domy,&isom,&qalfa,&qbeta);
//      // printf("range error = %f\n",range_error);
//      // printf("cond error = %f\n",max_adaptive_error * (long)(size *size));
//      if (s_size == min_size || range_error <= max_adaptive_error * (long)(s_size *s_size)) {
      
//         if(s_size > min_size)
//             pack(1,(long)0,fp);

//         if(abs(qalfa - zeroalfa) <= zero_threshold) {
//             qbeta = best_beta(atx,aty,s_size,0.0);
//             qalfa = zeroalfa;
//         }
//         pack(N_BITALFA, (long)qalfa, fp);
//         pack(N_BITBETA, (long)qbeta, fp);
//         if(isColor){
//           double um = u_mean(s_size,s_size,atx,aty);
//           double vm = v_mean(s_size,s_size,atx,aty);
//           pack(8,(long)um,fp);
//           pack(8,(long)vm,fp);
//         }
//         // printf("qalfa = %d\n",qalfa);
//         zero_alfa_transform ++;
//         if(qalfa != zeroalfa) {
//            zero_alfa_transform --;
//            pack(3, (long)isom, fp);
//            pack(bits_per_coordinate_h,(long)(domx / SHIFT),fp);
//            pack(bits_per_coordinate_w,(long)(domy / SHIFT),fp);
//         }

//         transforms ++;
//         coded += s_size * s_size ;
//         newper = (int) (((double) coded / (image_width * image_height)) * 100);
//         if(newper > oldper) {
//            bytes = pack(-2,(long)0,fp);
//            compress =(double) coded / (double)bytes;
//            bpp      = 8.0 /compress;
//            printf(" %d %% Coded, Compression %f:1, bpp %f\r", 
//                                       newper,compress,bpp);
//            oldper = newper;
//            fflush(stdout);
//         }
        
//      } else {

//         pack(1,(long)1,fp);

//         int h,hmax,vmax1,vmax2,i,j,y1,y2,Rp,x1,x2,x3,x4,Cp1,Cp2,v1,v2;
//         hmax=vmax1=vmax2=0;

//         Cp1=Cp2=atx;
//         Rp=aty;
//         for(j=aty;j<aty+s_size-1;j++)
//         {
//           y1=y2=0;
//           for(i=atx;i<atx+s_size;i++)
//           {
//             y1=y1+ range[j][i];
//             y2=y2+ range[j+1][i];
//           }
//           h=(min_2(j-aty+1,s_size-j+aty-1)) *abs(y1-y2);
//           if(h>hmax)
//           {
//             hmax=h;         
//             Rp=j;
//           }
//         }
//         for(i=atx;i<atx+s_size-1;i++)
//         {
//           x1=x2=x3=x4=0;
//           for(j=aty;j<Rp;j++)
//           {
//                 x1=x1+range[j][i];
//                 x2=x2+range[j][i+1];
//           }
//           v1=(min_2(i-atx+1,s_size-i+atx-1))*abs(x1-x2);
//           if(v1>vmax1)
//           {
//                 vmax1=v1;
//                 Cp1=i;
//           }
//           for(j=Rp;j<aty+s_size;j++)
//           {
//                 x3=x3+range[j][i];
//                 x4=x4+range[j][i+1];
//           }
//           v2=(min_2(i-atx+1 , s_size-i+atx-1))*abs(x3-x4);
//           if(v2>vmax2)
//           {
//                 vmax2=v2;
//                 Cp2=i;
//           }
//         }
//         // printf("%s\n","in adaptive part");
//         pack(3,(long)Rp-aty+1,fp);
//         pack(3,(long)Cp1-atx+1,fp);
//         pack(3,(long)Cp2-atx+1,fp);

//         adaptive_quadtree(atx,    aty,    min_2(Cp1-atx+1,Rp -aty+1));
//         adaptive_quadtree(Cp1+1,  aty,    min_2(atx+s_size-Cp1-1,Rp-aty+1));
//         adaptive_quadtree(atx,    Rp+1,   min_2(Cp2-atx+1,aty+s_size-Rp-1));
//         adaptive_quadtree(Cp2+1,  Rp+1,   min_2(atx+s_size-Cp2-1,aty+s_size-Rp-1));

        
//       }
//   }

//   if(x_size > s_size )
//      adaptive_quadtree(atx + s_size, aty, x_size - s_size, y_size);

//   if (y_size > s_size)
//     adaptive_quadtree(atx, aty + s_size, s_size, y_size - s_size);
// }