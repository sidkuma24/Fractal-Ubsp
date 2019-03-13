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
#include <string.h>

#include "nn_search.h"
#include "def.h"

#ifndef EXTERN
# define EXTERN extern
# define INIT(x)
   extern int ordering[][6];
   extern int mapping[][8];
#else
# define INIT(x) x

int ordering[][6] =  { {3,2,1,0,  3,0},
                       {2,3,1,0,  5,1},
                       {2,1,3,0,  5,2},
                       {3,1,2,0,  7,0},
                       {1,3,2,0,  1,1},
                       {1,2,3,0,  1,2},
                       {3,2,0,1,  3,1},
                       {2,3,0,1,  5,0},
                       {2,1,0,3,  2,2},
                       {3,1,0,2,  7,1},
                       {1,3,0,2,  1,0},
                       {1,2,0,3,  4,2},
                       {3,0,2,1,  3,2},
                       {2,0,3,1,  2,0},
                       {2,0,1,3,  2,1},
                       {3,0,1,2,  7,2},
                       {1,0,2,3,  4,1},
                       {1,0,3,2,  4,0},
                       {0,3,2,1,  6,2},
                       {0,2,3,1,  6,1},
                       {0,2,1,3,  6,0},
                       {0,3,1,2,  0,2},
                       {0,1,3,2,  0,1},
                       {0,1,2,3,  0,0}
                    };

int mapping[][8]  =  { {0, 1, 2, 3, 4, 5, 6, 7},
                       {2, 0, 3, 1, 7, 6, 4, 5},
                       {1, 3, 0, 2, 6, 7, 5, 4},
                       {3, 2, 1, 0, 5, 4, 7, 6},
                       {4, 7, 6, 5, 0, 3, 2, 1},
                       {5, 6, 7, 4, 3, 0, 1, 2},
                       {6, 4, 5, 7, 1, 2, 0, 3},
                       {7, 5, 4, 6, 2, 1, 3, 0}
                    };
#endif

typedef unsigned char PIXEL;


struct code_book {
            short ptr_x;
            short ptr_y;
            short isom;
            double sum;
            double sum2;
            double sum3;
            double sum4;
         };

struct c {
	   short ptr_x;
	   short ptr_y;
	   double sum;
	   double sum2;
     double sum3;
     double sum4;
     double var;
     int iso;

	   struct c *next;
	 };


struct t_node {
           double rx;
           double ry;
           double rrx;
           double rry;
           short size;
           int  x_size;
           int y_size;
           double dx;
           double dy;
           short um;  // for color image: uchannel mean
           short vm;  // for color image: vchannel mean
           short sym_op;
           double alfa,alfa1,alfa2, beta, rmean;        
           int qdx,qdy,qalfa1,qalfa2,qalfa, qbeta,qrmean;
           struct t_node *next;
         } ;

EXTERN struct t_node *trans,fractal_code;

EXTERN int min;
EXTERN int max;
EXTERN int lev;

EXTERN FILE *input;
EXTERN FILE *output;
EXTERN int iterations   INIT(= 10);
EXTERN int postproc     INIT(= 0);
EXTERN int quality      INIT(= 0);
EXTERN int piramidal    INIT(= 1);
EXTERN int iterDec2     INIT(= 0);
EXTERN int raw_format   INIT(= 0);
EXTERN int display      INIT(= 0);
EXTERN double zoom      INIT(= 1.0);

EXTERN PIXEL **image1;
EXTERN PIXEL **init_image;
EXTERN double max_std_arr[8];
EXTERN double max_ent_arr[8];
EXTERN int final_max_std INIT(=0);
EXTERN double final_max_ent INIT(=0);
EXTERN int final_max_ent_q INIT(=0);
EXTERN int save_count INIT(=0);

EXTERN struct c ***class_polar[8];
EXTERN struct c *class_fisher[8][3][24];
// EXTERN struct c **class_std[8][3];
EXTERN struct c *class_std[8][3][255];
EXTERN struct c **class_entropy[8][3];
EXTERN struct c *class_basicFIC[8][3];
// EXTERN struct c *class_fisher[8][3][240];
EXTERN struct c *class_hurtgen[8][16][24];

EXTERN kdtree  ***class_polar_saupe[8];
EXTERN int     **item_per_class[8];
EXTERN struct code_book  ***c_book[8];
EXTERN float ****f_vect[8];  

EXTERN struct code_book  *codebook[8];
EXTERN float **f_vectors[8];  
EXTERN kdtree *kd_tree[8];
EXTERN int feat_vect_dim[8];
EXTERN int average_factor[8];
EXTERN int n_features INIT(= 16);
EXTERN int shrunk_factor_mc INIT(= 1);
EXTERN int shrunk_factor_saupe INIT(= 0);

EXTERN int N_BITALFA  INIT(=  5);
EXTERN int N_BITBETA  INIT(=  7);

/* for nonlinear error checking */
EXTERN int N_BITALFA1 INIT(=  5);
EXTERN int N_BITALFA2  INIT(=  5);
EXTERN int N_BITBETA2 INIT(=  8);
EXTERN int N_BITRMEAN INIT(=  7);


EXTERN double MAX_ALFA  INIT(=  1.0);
EXTERN int MAX_BITS  INIT(=  3);
EXTERN double MAX_ALFA1  INIT(=  1.5);
EXTERN double MAX_ALFA2 INIT(=  0.015);
EXTERN int min_size INIT(= 4);
EXTERN int max_size INIT(= 16);
EXTERN int image_width     INIT(= 512);
EXTERN int image_height    INIT(= 512);
EXTERN int virtual_size    INIT(= 512);
EXTERN int SHIFT    INIT(= 4);

EXTERN double T_ENT  INIT(= 8.0);
EXTERN double T_VAR  INIT(= 1000000.0);
EXTERN double T_RMS  INIT(= 8.0);
EXTERN double adapt INIT(= 1.0);

EXTERN double eps INIT(= 2.0);
EXTERN int matches INIT(= 50);
EXTERN int n_p_class INIT(= 50);
EXTERN int livelli  INIT(= 3);
EXTERN int zero_threshold INIT(= 0);
EXTERN int best     INIT(= 0);
EXTERN int full_first_class INIT(= 0);
EXTERN int full_second_class INIT(= 0);
EXTERN int qtree    INIT(= 0);
EXTERN FILE *fp;

EXTERN int bits_per_coordinate_w;
EXTERN int bits_per_coordinate_h;

EXTERN int zero_alfa_transform INIT(= 0);
EXTERN long comparisons INIT(= 0);
EXTERN int transforms INIT(= 0);
EXTERN int zeroalfa;

EXTERN char filein[100];
EXTERN char fileout[100];

EXTERN PIXEL **image;
EXTERN PIXEL **image_uch;
EXTERN PIXEL **image_vch;
EXTERN PIXEL **qtt;
EXTERN PIXEL **hv;
EXTERN double **contract;
EXTERN double **range;
EXTERN double **range_tmp;
EXTERN double **flip_range;
EXTERN double (*Coding) (int, int, int, int*, int*, int*, int*, int*);
EXTERN double (*testing_Coding) (int, int, int, int*, int*, int*, int*, int*);
EXTERN double (*LumInv_Coding) (int, int, int, int*, int*, int*, int*, int*);
EXTERN double (*Nonlinear_Coding) (int, int, int, int*, int*, int*, int*, int*,int* );
EXTERN void (*Indexing) (int, int);
EXTERN int method INIT(= MassCenter);
EXTERN int isColor INIT(= 0);
EXTERN int isNonlinear INIT(=0);
EXTERN int isTesting INIT(=0);
EXTERN int isLumInv INIT(=0);
EXTERN int isCovar2 INIT(=0);
EXTERN int zero_alfa_count INIT(=0);
EXTERN double BPP INIT(=0);
EXTERN double COMPRESS INIT(=0);
EXTERN double max_adaptive_error INIT(=0.1);
EXTERN double range_error INIT(=0);

/* for adaptive quadtree */
extern int DOMAIN_X_BITS;
extern int DOMAIN_Y_BITS;
extern int MIN_ADAP_D_BITS INIT(=  0); // min size = 2;
extern int MAX_ADAP_D_BITS INIT(=  3); // max size = 64; i.e 2^6
EXTERN struct c* adaptive_fisher_class[64][64][24];
EXTERN void (*adaptiveIndexing) (int, int);
extern int MIN_ADAP_R_BITS INIT(=  MIN_ADAP_R_BITS);
extern int MAX_ADAP_R_BITS INIT(=  MAX_ADAP_D_BITS );

EXTERN double (*adaptiveCoding) (int, int, int, int,int*, int*, int*, int*, int*);

EXTERN int **cum_range;
EXTERN int **cum_range2;
EXTERN int **L1_var_limits;
EXTERN int *L1_class_width;
EXTERN int N_L1_CLASSES INIT(=24); 
EXTERN int N_L2_CLASSES INIT(=24); 

EXTERN int partition_type INIT(= Quadtree);

EXTERN double max_error2 INIT(=100);


/* HV parition */
EXTERN void (*HV_Indexing) (int, int);
EXTERN struct c* HV_fisher_class[64][64][24];
EXTERN double (*HV_Coding) (int, int, int, int,int*, int*, int*, int*, int*);
EXTERN int isHV INIT(=0);


/* progressive decoding */
EXTERN int N_BITS INIT(=7);

extern int class_count[24];
