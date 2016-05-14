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
#include "prot.h"


void ComputeMc(double **block,int size,double *x,double *y,int s)
{
  double row = 0.0;
  double col = 0.0;
  double a = 0.0;
  double b = 0.0;
  double mass  = 0.0;
  register int i,j;
  
  for(i=0;i<size;i++) {
    a = 0.0;
    b = 0.0;
    for(j=0;j<size;j++) {
      a += block[i][j] ;
      b += block[j][i] ;
      mass += block[j][i];
    }
    row += a * (i+1);
    col += b * (i+1); 
  }
  if(mass != 0.0) {
      *y = (size+1)/2.0 - row / mass ; 
      *x = col / mass - (size+1)/2.0 ; 

  }
  else {
   *x = 0.0;
   *y = 0.0;
  }
}


int match(int *iso_1)
{
  int i,j,flag;

  for(j=0; j< 24; j++) {
    flag =1;
    for(i=0; i< 4; i++)
      if(iso_1[i] != ordering[j][i]) {
        flag = 0;
        break;
      }
    if (flag == 1) return j;
  }
  fatal("\n Error");
  return(-1);
}


int variance_class(int size,double **block)
{
  int i,j;
  double a[4] = {0.0,0.0,0.0,0.0};
  int order[4];

  a[0] =  variance_2(size/2, block, 0, 0);
  a[1] =  variance_2(size/2, block, 0, size/2);
  a[2] =  variance_2(size/2, block, size/2, 0);
  a[3] =  variance_2(size/2, block, size/2, size/2);

  for (i=0; i<4; ++i)  order[i] = i ;

  for (i=2; i>=0; --i)
  for (j=0; j<=i; ++j)
     if (a[j]<a[j+1]) {
       swap(order[j], order[j+1],int)
       swap(a[j], a[j+1],double)
     }

  return  match(order);
}


int hurtgen_class(int size,double **block)
{
  int i,j;
  int clas = 0;
  double quadrant_mean[4] = {0.0,0.0,0.0,0.0};
  double mean;

  for(i=0;i<size/2;i++)
  for(j=0;j<size/2;j++) {
     quadrant_mean[0] += block[i][j];
     quadrant_mean[1] += block[i][j+size/2];
     quadrant_mean[2] += block[i+size/2][j];
     quadrant_mean[3] += block[i+size/2][j+size/2];
  }

  mean = (quadrant_mean[0] + quadrant_mean[1] + 
	           quadrant_mean[2] + quadrant_mean[3]) / (size * size);

  quadrant_mean[0] /= (size/2 * size/2);
  quadrant_mean[1] /= (size/2 * size/2);
  quadrant_mean[2] /= (size/2 * size/2);
  quadrant_mean[3] /= (size/2 * size/2);

  for(i=0;i<4;i++) {
    clas <<=1;
    if(mean < quadrant_mean[i])
       clas |= 1;
  }
  return clas;
}


void newclass(int size,double **block,int *isom,int *clas)
{
  int i,j;
  int index;
  double a[4] = {0.0,0.0,0.0,0.0};
  int order[4];


  for(i=0;i<size/2;i++)
  for(j=0;j<size/2;j++) {
     a[0] += block[i][j];
     a[1] += block[i][j+size/2];
     a[2] += block[i+size/2][j];
     a[3] += block[i+size/2][j+size/2];
  }

  for (i=0; i<4; ++i)  order[i] = i
        ;

  for (i=2; i>=0; --i)
  for (j=0; j<=i; ++j)
     if (a[j]<a[j+1]) {
       swap(order[j], order[j+1],int)
       swap(a[j], a[j+1],double)
     }

  index = match(order);

  *isom = ordering[index][4];
  *clas = ordering[index][5];
}


void contraction(double **t,PIXEL **fun,int s_x,int s_y)
{
  double tmp;
  int i,j,k,w,s,z;

  for(i=s_x;i< image_height - s_x;i+=2)
  for(j=s_y;j< image_width  - s_y;j+=2) {
     tmp=0.0;
     for(k=0; k < 2; k++) {
       s = 0;
       if((i+ k) >= image_height) 
	  s = 1;
       for(w=0; w < 2; w++) {
         z = 0;
         if((j+ w) >= image_height) 
           z = 1;
         tmp +=(double) fun[i+k-s][j+w-z];
         t[i >> 1][j >> 1] = tmp/4.0;
       }
     }
  }
}


int ShrunkBlock(double **block,int size, int srunk_fact)
{
  double tmp;
  int i,j,k,w;
  int shift;

  shift = (int) rint(log((double) srunk_fact) / log(2.0));

  for(i=0;i< size ;i+=srunk_fact)
  for(j=0;j< size ;j+=srunk_fact) {
     tmp=0;
     for(k=0; k < srunk_fact; k++)
     for(w=0; w < srunk_fact; w++)
        tmp += block[i+k][j+w];
     block[i >> shift][j >> shift] = tmp / (double) (srunk_fact*srunk_fact);
  }
  return(size / srunk_fact);
}


void ComputeSaupeVectors(double **block, int size,int index, float *vector)
{
  double sum = 0.0;
  double var = 0.0;
  double s, v ;
  register int j,i;
  register int k=0;
  int newsize;

  if(average_factor[index] > 1)
     newsize = ShrunkBlock(block,size,average_factor[index]);
  else
     newsize = size;

  for(i=0;i<newsize;i++)
  for(j=0;j<newsize;j++) {
     sum += block[i][j];
     var += block[i][j] * block[i][j];
  }

  s = sum / (double)(newsize*newsize);
  v = sqrt(var - sum*sum / (double)(newsize*newsize));

  for(i=0;i<newsize;i++)
  for(j=0;j<newsize;j++,k++)
    vector[k] = (block[i][j] - s ) / v;

}


void ComputeMcVectors(double **block,double **block_tmp,
                                    int size, int index,double *vector)
{
  double mean;
  double xx,yy;
  double theta;
  double sum = 0.0;
  double **block_1;
  double **block_2;
  register int x,y,i;
  int newsize;

  block_1 = block;
  block_2 = block_tmp;

  if(average_factor[index] > 1)
     newsize = ShrunkBlock(block,size,average_factor[index]);
  else
     newsize = size;
 
  ComputeMc(block_1,newsize, &xx, &yy, 0); /* Compute the MC of the block     */
  theta = atan2(yy, xx);                   /* and put it in polar coordinates */
  if(theta >= 0.0) 
     vector[0] = theta;
  else  
     vector[0] = TWOPI + theta;

  for(x=0;x< newsize;x++)
  for(y=0;y< newsize;y++) {
     sum += block_1[x][y];
  }

  for(i=1;i<2;i++) {                    /* Compute other features by       */
     mean = sum / (newsize * newsize);  /* transforming the original block */
     sum = 0.0;
     for(x=0;x< newsize;x++)
     for(y=0;y< newsize;y++) {
        block_2[x][y] = (block_1[x][y] - mean) * (block_1[x][y] - mean);
        sum += block_2[x][y];
     }
     ComputeMc(block_2,newsize, &xx, &yy,i);
     theta = atan2(yy, xx);     
     if(theta >= 0.0) 
        vector[i] = theta;
     else  
        vector[i] = TWOPI + theta;

     swap(block_1, block_2, double **)
  }

}


void Saupe_FisherIndexing(int size,int s)
{
  int i,j,k;
  int count = 0;
  int cbook_size;
  int dim_vectors;
  int clas,iso;
  double sum, sum2;
  double **dom_tmp,**domi, **flip_domi;
  register double pixel;
  register int x,y;

  cbook_size = (1+image_width / SHIFT) * (1+image_height / SHIFT);

  dim_vectors = feat_vect_dim[(int) rint(log((double)(size))/log (2.0))];
  matrix_allocate(f_vectors[s],dim_vectors,cbook_size,float)
  codebook[s]=(struct  code_book*) malloc(cbook_size*sizeof(struct code_book));

  matrix_allocate(domi,size,size,double)
  matrix_allocate(flip_domi,size,size,double)
  matrix_allocate(dom_tmp,size,size,double)

  for(i=0;i< image_height - 2 * size +1 ;i+= SHIFT) {
    for(j=0;j< image_width  - 2 * size +1 ;j+= SHIFT,count ++) {
      k=0;
      sum  = 0.0;
      sum2 = 0.0;
      for(x=0;x< size;x++)
      for(y=0;y< size;y++) {
        pixel = contract[x+(i>>1)][y+(j>>1)];
        sum  += pixel;
        sum2 += pixel * pixel;
        domi[x][y] = pixel;
      }
                                      /* Compute the symmetry operation which */
      newclass(size,domi,&iso,&clas); /* brings the domain in the canonical   */
                                      /* orientation                          */
                                      
      flips(size,domi,flip_domi,iso); 

      ComputeSaupeVectors(flip_domi,size,s,f_vectors[s][count]);

      codebook[s][count].sum   = sum;
      codebook[s][count].sum2  = sum2;
      codebook[s][count].ptr_x = i;
      codebook[s][count].ptr_y = j;
      codebook[s][count].isom  = iso;
    }
    printf(" Extracting [%d] features (Saupe)  domain (%dx%d)  %d \r",
	                         	    dim_vectors, size,size,count) ;
    fflush(stdout);

  }
  printf("\n Building Kd-tree... ") ;
  fflush(stdout);
  kd_tree[s] = kdtree_build(f_vectors[s],count-1,dim_vectors);
  printf("Done\n") ;
  fflush(stdout);

  free(domi[0]);
  free(dom_tmp[0]);
  free(flip_domi[0]);
}



void MassCenterIndexing(int size,int s)
{
  int i,j,k;
  int count = 0;
  int cx,cy;
  struct c *node; 
  double sum,sum2;

  double vector[16];

  double **dom_tmp,**domi,**flip_domi;
  register double pixel;
  register int x,y;

  matrix_allocate(domi,size,size,double)
  matrix_allocate(dom_tmp,size,size,double)
  matrix_allocate(flip_domi,size,size,double)

  for(i=0;i< image_height - 2 * size +1 ;i+= SHIFT) {
    for(j=0;j< image_width  - 2 * size +1 ;j+= SHIFT,count ++) {
      k=0;
      sum  = 0.0;
      sum2 = 0.0;
      for(x=0;x< size;x++)
      for(y=0;y< size;y++) {
        pixel = contract[x+(i>>1)][y+(j>>1)];
        sum  += pixel;
        sum2 += pixel * pixel;
        domi[x][y] = pixel;
      }

      ComputeMcVectors(domi,dom_tmp,size,s,vector);

      cx = (int) (vector[0] / TWOPI * n_p_class ); 
      cy = (int) (vector[1] / TWOPI * n_p_class ); 

      node = (struct c *) malloc(sizeof(struct c));
      node -> ptr_x = i;
      node -> ptr_y = j;
      node -> sum   = sum;
      node -> sum2  = sum2;
      node -> next  = class_polar[s][cx][cy];
      class_polar[s][cx][cy] = node;
    }
    printf(" Extracting features (Mass Center) domain (%dx%d)  %d \r",
                                                              size,size,count) ;
    fflush(stdout);

  }

  printf("\n") ;

  free(domi[0]);
  free(dom_tmp[0]);
  free(flip_domi[0]);
}


void Mc_SaupeIndexing(int size,int s)
{
  int i,j,k,h;
  int count = 0;
  int cx,cy;
  struct c *node; 
  int dim_vectors;
  struct c *pointer;
  double sum,sum2;
  double vector[16];
  double **dom_tmp,**domi,**flip_domi;
  register double pixel;
  register int x,y;

  matrix_allocate(domi,size,size,double)
  matrix_allocate(dom_tmp,size,size,double)
  matrix_allocate(flip_domi,size,size,double)

  for(i=0;i< image_height - 2 * size +1 ;i+= SHIFT) {
    for(j=0;j< image_width - 2 * size +1 ;j+= SHIFT,count ++) {
      k=0;
      sum  = 0.0;
      sum2 = 0.0;
      for(x=0;x< size;x++)
      for(y=0;y< size;y++) {
        pixel = contract[x+(i>>1)][y+(j>>1)];
        sum  += pixel;
        sum2 += pixel * pixel;
        domi[x][y] = pixel;
      }

      ComputeMcVectors(domi,dom_tmp,size,s,vector);

      cx = (int) (vector[0] / TWOPI * n_p_class ); 
      cy = (int) (vector[1] / TWOPI * n_p_class ); 

      node = (struct c *) malloc(sizeof(struct c));
      node -> ptr_x = i;
      node -> ptr_y = j;
      node -> sum   = sum;
      node -> sum2  = sum2;
      node -> next  = class_polar[s][cx][cy];
      class_polar[s][cx][cy] = node;

      item_per_class[s][cx][cy] ++;
    }
    printf(" Extracting features (Mass Center) domain (%dx%d)  %d \r",
                                                            size,size,count) ;
    fflush(stdout);

  }
  printf("\n") ;

  dim_vectors = feat_vect_dim[(int) rint(log((double) size)/log (2.0))];

  printf(" Extracting features (Saupe) domain (%dx%d) building Kd-trees",
                                                                size,size) ;
  fflush(stdout);

  for(h=0;h<n_p_class;h++)
  for(i=0;i<n_p_class;i++) {
    c_book[s][h][i] =(struct code_book*) malloc(item_per_class[s][h][i] *
                                                       sizeof(struct code_book));
    matrix_allocate(f_vect[s][h][i],dim_vectors,item_per_class[s][h][i],float)
    pointer = class_polar[s][h][i];
    count = 0;
    while(pointer != NULL) {
      for(x=0;x< size;x++)
      for(y=0;y< size;y++) 
        domi[x][y] = contract[x+(pointer->ptr_x>>1)][y+(pointer->ptr_y>>1)];

      ComputeSaupeVectors(domi,size,s,f_vect[s][h][i][count]);

      c_book[s][h][i][count].sum   = pointer -> sum;
      c_book[s][h][i][count].sum2  = pointer -> sum2;
      c_book[s][h][i][count].ptr_x = pointer -> ptr_x;
      c_book[s][h][i][count].ptr_y = pointer -> ptr_y;

      count ++;
      pointer = pointer -> next;
    }

    if(count > 0)                                             
      class_polar_saupe[s][h][i] = kdtree_build(f_vect[s][h][i],count,dim_vectors);
  }
  printf("\n") ;
  free(domi[0]);
  free(flip_domi[0]);
  free(dom_tmp[0]);
}



void SaupeIndexing(int size,int s)
{
  int i,j,k;
  int count = 0;
  int cbook_size;
  int dim_vectors;
  double sum,sum2;
  double **dom_tmp,**domi;
  register double pixel;
  register int x,y;

  cbook_size = (1+image_width / SHIFT) * (1+image_height / SHIFT);
  dim_vectors = feat_vect_dim[(int) rint(log((double) size)/log (2.0))];

  matrix_allocate(f_vectors[s],dim_vectors,cbook_size,float)
  codebook[s]=(struct  code_book*) malloc(cbook_size*sizeof(struct code_book));
  matrix_allocate(domi,size,size,double)
  matrix_allocate(dom_tmp,size,size,double)

  for(i=0;i< image_height - 2 * size +1 ;i+= SHIFT) {
    for(j=0;j< image_width - 2 * size +1 ;j+= SHIFT,count ++) {
      k=0;
      sum  = 0.0;
      sum2 = 0.0;
      for(x=0;x< size;x++)
      for(y=0;y< size;y++) {
        pixel = contract[x+(i>>1)][y+(j>>1)];
        sum  += pixel;
        sum2 += pixel * pixel;
        domi[x][y] = pixel;
      }

      ComputeSaupeVectors(domi,size,s,f_vectors[s][count]);

      codebook[s][count].sum  = sum;
      codebook[s][count].sum2 = sum2;
      codebook[s][count].ptr_x = i;
      codebook[s][count].ptr_y = j;
    }
    printf(" Extracting [%d] features (Saupe) domain (%dx%d)  %d \r",
		                               dim_vectors, size,size,count) ;
    fflush(stdout);

  }
  printf("\n Building Kd-tree... ") ;
  fflush(stdout);
  kd_tree[s]= kdtree_build(f_vectors[s],count-1,dim_vectors);
  printf("Done\n") ;
  fflush(stdout);

  free(domi[0]);
  free(dom_tmp[0]);

}


void HurtgenIndexing(int size,int s)
{
  int i,j,k,h;
  int count = 0;
  int m_class_1, m_class_2;
  double sum,sum2;
  double **domi;
  register double pixel;
  register int x,y;
  struct c *node;

  matrix_allocate(domi,size,size,double)

  for(i=0;i< image_height - 2 * size +1 ;i+= SHIFT) {
    for(j=0;j< image_width - 2 * size +1 ;j+= SHIFT) {
      count ++;
      k=0;
      sum  = 0.0;
      sum2 = 0.0;
      for(x=0;x< size;x++)
      for(y=0;y< size;y++) {
        pixel = contract[x+(i>>1)][y+(j>>1)];
        sum  += pixel;
        sum2 += pixel * pixel;
        domi[x][y] = pixel;
      }

      m_class_1 = hurtgen_class(size,domi);
      m_class_2 = variance_class(size,domi);

      node = (struct c *) malloc(sizeof(struct c));
      node -> ptr_x = i;
      node -> ptr_y = j;
      node -> sum   = sum;
      node -> sum2  = sum2;
      node -> next  = class_hurtgen[s][m_class_1][m_class_2];
      class_hurtgen[s][m_class_1][m_class_2] = node;
    }
    printf(" Classification (Hurtgen) domain (%dx%d)  %d \r",size,size,count) ;
    fflush(stdout);
  }

  for(i=0;i<16;i++)     /* Find a not empty class */
  for(j=0;j<24;j++)
     if(class_hurtgen[s][i][j] != NULL)
          goto out_loops;


  out_loops:

  for(k=0;k<16;k++)   /* Make sure no class is empty */
  for(h=0;h<24;h++)
     if(class_hurtgen[s][k][h] == NULL)
       class_hurtgen[s][k][h] = class_hurtgen[s][i][j];


  printf("\n");
  free(domi[0]);

}



void FisherIndexing(int size,int s)
{
  int i,j,k,h;
  int count = 0;
  int iso, clas, var_class;
  double sum,sum2;
  double **domi,**flip_domi;
  register double pixel;
  register int x,y;
  struct c *node;

  matrix_allocate(domi,size,size,double)
  matrix_allocate(flip_domi,size,size,double)

  for(i=0;i< image_height - 2 * size +1 ;i+= SHIFT) {
    for(j=0;j< image_width - 2 * size +1 ;j+= SHIFT ) {
      count ++;
      k=0;
      sum  = 0.0;
      sum2 = 0.0;
      for(x=0;x< size;x++)
      for(y=0;y< size;y++) {
        pixel = contract[x+(i>>1)][y+(j>>1)];
        sum  += pixel;
        sum2 += pixel * pixel;
        domi[x][y] = pixel;
      }
                                      /* Compute the symmetry operation which */
      newclass(size,domi,&iso,&clas); /* brings the domain in the canonical   */
                                      /* orientation                          */
      flips(size,domi,flip_domi,iso); 

      var_class = variance_class(size,flip_domi);

      node = (struct c *) malloc(sizeof(struct c));
      node -> ptr_x = i;
      node -> ptr_y = j;
      node -> sum   = sum;
      node -> sum2  = sum2;
      node -> iso  = iso;
      node -> next  = class_fisher[s][clas][var_class];
      class_fisher[s][clas][var_class] = node;
    }
    printf(" Classification (Fisher) domain (%dx%d)  %d \r",size,size,count) ;
    fflush(stdout);

  }

  for(i=0;i<3;i++)    /* Find a not empty class */
  for(j=0;j<24;j++)
    if(class_fisher[s][i][j] != NULL)
         goto out_loops;


  out_loops:

  for(k=0;k<3;k++)  /* Make sure no class is empty */
  for(h=0;h<24;h++)
    if(class_fisher[s][k][h] == NULL)
      class_fisher[s][k][h] = class_fisher[s][i][j];


  printf("\n");

  free(domi[0]);
  free(flip_domi[0]);
}


void ComputeFeatVectDimSaupe()
{
  int i;
  int size;

  for(i=1;i<=(int) rint(log((double) max_size)/ log(2.0)); i++) {
    size = (int) rint(pow(2.0,(double)i));
    if(shrunk_factor_saupe) {
      n_features = (size*size) / (shrunk_factor_saupe * shrunk_factor_saupe);
      if(n_features < 4)
        n_features = 4;
    }
    if((size*size)  > n_features ) { 
       average_factor[i] = size / sqrt(n_features);
       feat_vect_dim[i] = n_features;
    }
    else {
       average_factor[i] = 1;
       feat_vect_dim[i] = size*size;
    }
  }
}


void ComputeAverageFactorMc()   
{                             
  int i;
  int size;

  for(i=1;i<=(int) rint(log((double) max_size)/ log(2.0)); i++) {
      size = (int) rint(pow(2.0,(double)i));
      if(shrunk_factor_mc) {
          average_factor[i] = shrunk_factor_mc;
          if((size*size)/(shrunk_factor_mc*shrunk_factor_mc) < 16) 
          average_factor[i] = 1;
       }
       else {
          average_factor[i] = 1;
       }
   }
}



