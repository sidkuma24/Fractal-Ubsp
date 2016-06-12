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


int main(int argc, char **argv)
{
  int int_max_alfa;

  printf("%s",triangle_4);

  getopt_dec(argc, argv);

  if ((input = fopen(filein, "r")) == NULL)
    return -1;

  unpack(-2,input);    /*Initialize unpack */

  N_BITALFA      = (int)unpack(4,input);
  N_BITBETA      = (int)unpack(4,input);
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
  read_transformations(0,0,virtual_size);
  printf("done\n");
  fflush(stdout);

    return 0;
}

void read_transformations(int atx,int aty,int size)
{ 
  int qalfa,qbeta;
  double alfa,beta, um,vm;


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
          trans->dx = SHIFT * (int)unpack(bits_per_coordinate_h,input);
          trans->dy = SHIFT * (int)unpack(bits_per_coordinate_w,input);
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

