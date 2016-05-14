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

#include "def.h"
#include "globals.h"
#include "nn_search.h"
#include "prot.h"



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


double FisherCoding(int atx,int aty,int size,int *xd,int *yd,int *is,
                                                       int* qalf,int *qbet)
{
  int x,y,qalfa ,qbeta,i,j;
  int tip,counter,isometry ;
  int isom,clas, var_class;
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


