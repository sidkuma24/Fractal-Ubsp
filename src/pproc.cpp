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

// FILE * output;
int globalcounter = 0;

int main(int argc, char **argv)
{
  int int_max_alfa;

  printf("%s",triangle_4);

  getopt_dec(argc, argv);

  if ((input = fopen(filein, "r")) == NULL){
    printf("error: Cannot open file %s\n", filein);
    return -1;
  }

  // Open the output file
  int i = 0;
  char fileout[100];
  while (filein[i] != '.') {
      fileout[i] = filein[i];
      i++;
  }
  strcat(fileout, ".ifsp");

  if ((output = fopen(fileout, "w")) == NULL){
    printf("error: Cannot open file %s\n", fileout);
    return -1;
  }
    
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
  printf("Minsize: %d\n", min_size);
  
  /* Write the header to the output file */
  pack(4, N_BITALFA    , output);
  pack(4, N_BITBETA    , output);
  pack(7, min_size     , output);
  pack(7, max_size     , output);
  pack(6, SHIFT        , output);
  pack(12, image_width , output);
  pack(12, image_height, output);
  pack(8, int_max_alfa , output);
  pack(1, isColor      , output);


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
  pack(32, 420, output);

  read_transformations(0,0,virtual_size);
  printf("done\n");
  fflush(stdout);

  // int bit_depth = 0;
  // while(bit_depth < 8){       /* TODO: 8 is for color bits, change this to the maximum bit length parameter */
  //   printf("%s\n", "Press [Enter] to write the MSB of parameters:");
  //   getchar();
  //   write_details(bit_depth);
  //   bit_depth++;
  // }

  printf("No of transforms: %d\n", globalcounter);
  // Close the files
  pack(-1, (long)0, output);
  pack(-2, (long)0, output);
  fclose(output);
  fclose(input);
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
      printf("1\n");
      pack(1, (long)1, output);
      read_transformations(atx,aty,size/2);
      read_transformations(atx+size/2,aty,size/2);
      read_transformations(atx,aty+size/2,size/2);
      read_transformations(atx+size/2,aty+size/2,size/2);
  } else {
      /* Read the trasformation */ 
      //printf("inside else\n"); 
      // globalcounter++;
      printf("0\n");
      pack(1, (long)0, output);
      // trans->next = (struct t_node *) malloc(sizeof(struct t_node ));
      // trans       = trans->next; 
      // trans->next = NULL;
      // trans->qalfa       = (int)unpack(N_BITALFA,  input);
      // trans->qbeta       = (int)unpack(N_BITBETA,  input);
      // if(isColor){
      //   trans->um = (int)unpack(8,input);
      //   trans->vm = (int)unpack(8,input);
      // }

      // if(trans->qalfa != zeroalfa) {      
      //     trans-> sym_op = (int)unpack(3, input);
      //     trans->qdx = (int)unpack(bits_per_coordinate_h,input);
      //     trans->qdy = (int)unpack(bits_per_coordinate_w,input);
      // } else {
      //     trans->sym_op = 0;
      //     trans->qdx  = 0;
      //     trans->qdy = 0;
      // }
      // pack(3, (long)trans->sym_op, output);
      // pack(bits_per_coordinate_h, (long)trans->qdx, output);
      // pack(bits_per_coordinate_w, (long)trans->qdy, output);
      // printf("%d %d\n", trans->qdx, trans->qdy);

  }
}

void write_details(int bit_depth){
  trans = &fractal_code;
  int bits, value;
  while(trans->next){
    value = 0;
    bits = 0;
    trans = trans->next;
    if(bit_depth < N_BITALFA){
      value |= trans->qalfa & (1<<(N_BITALFA - bit_depth -1));
      value <<= 1;
      bits++;
    }
    if(bit_depth < N_BITBETA){
      value |= trans->qbeta & (1<<(N_BITBETA - bit_depth -1));
      value <<= 1;        
      bits++;
    }
    if(isColor){
      if (bit_depth < 8)
      {
        value |= trans->um & (1<<(8 - bit_depth -1));
        value <<= 1;        
        value |= trans->vm & (1<<(8 - bit_depth -1));
        value <<= 1;        
        bits += 2;
      }
    }
    pack(bits, (long)value, output);
    printf("%d\n", value);
  }
}