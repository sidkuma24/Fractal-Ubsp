/******************************************************************************
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



double HurtgenCoding(int,int,int,int *,int *,int *, int *,int *);
double SaupeCoding(int,int,int,int *,int *,int *, int *,int *);
double FisherCoding(int,int,int,int *,int *,int *, int *,int *);
double STDCoding(int,int,int,int *,int *,int *, int *,int *);
double COVCoding(int,int,int,int *,int *,int *, int *,int *);
double EntropyCoding(int,int,int,int *,int *,int *, int *,int *);
double AdaptiveSearch_FisherCoding(int,int,int,int *,int *,int *, int *,int *);
double CovClass_AdaptiveSearch_FisherCoding(int,int,int,int *,int *,int *, int *,int *);
double Nonlinear_FisherCoding(int,int,int,int *,int *,int *, int *,int *, int *);
double LumInv_FisherCoding(int,int,int,int *,int *,int *, int *,int *);
double testing_FisherCoding(int,int,int,int *,int *,int *, int *,int *);
double LumInv_FisherCoding2(int,int,int,int *,int *,int *, int *,int *);
double BasicFIC_Coding(int,int,int,int *,int *,int *, int *,int *);

double Mc_SaupeCoding(int,int,int,int *,int *,int *,int *,int *);
double MassCenterCoding(int,int,int,int *,int *,int *,int *,int *);
double Saupe_FisherCoding(int,int,int,int *,int *,int *,int *,int *);
double entropy(int, int,int,int);
double variance(int, int,int,int);
double variance_2(int, double **, int, int);
double entropy_2(int, double **, int, int);
double energy_coeff(int, double **, int, int);
int EnergyCoeff_class(int ,double **);
int std_class(int size, double **block);
int ent_class(int size, double **block);
void findMaxStd(int size, int s);
void findMaxEnt(int size, int s);

double u_mean(int, int,int,int);
double v_mean(int, int,int,int);
int hurtgen_class(int ,double **);
int pack(int , long , FILE *);
int variance_class(int ,double **);
kdtree *kdtree_build(float **,int,int);
void ComputeAverageFactorMc();
void ComputeFeatVectDimSaupe();
void ComputeMcVectors(double **,double **,int,int,double *);
void ComputeSaupeVectors(double **,int,int,float *);
void FisherIndexing(int ,int);
void BasicFIC_Indexing(int ,int);
void STDIndexing(int ,int);
void COVIndexing(int ,int);
void EntropyIndexing(int ,int);
void EnergyCoeff_FisherIndexing(int ,int);
void HurtgenIndexing(int ,int);
void MassCenterIndexing(int ,int);
void SaupeIndexing(int ,int);
void Mc_SaupeIndexing(int ,int);
void Saupe_FisherIndexing(int ,int);
void contraction(double **,PIXEL **,int,int);
void fatal(char *);
void flips(int ,double **,double **,int);
int kdtree_search(float *,float **,int,kdtree *,float,int,int*);
void ComputeMc(double **,int ,double *,double *, int);
void newclass(int,double **,int *,int *);
void getopt_enc(int,char **);
void getopt_dec(int,char **);
void quadtree(int ,int,int ,double ,double,double);
void quadtree_2(int ,int,int ,double ,double,double);
void testing_quadtree(int ,int,int ,double ,double,double);
void Nonlinear_quadtree(int ,int,int ,double ,double,double);
void LumInv_quadtree(int ,int,int ,double ,double,double);
void readimage_raw(char *);
void readimage_pgm(char *,int *,int *);
void help_enc();
int  bitlength(unsigned long val);
void read_transformations_2(int atx,int aty,int size);


long unpack(int size, FILE *fin);
void write_details(int bit_depth);
void read_details(int bit_depth);
void read_transformations(int atx,int aty,int size);
void read_transformations_new_init_image(int atx,int aty,int size);
void read_transformations_testing(int atx,int aty,int size);
void read_transformations_LumInv(int atx,int aty,int size);
void read_initial_transformations(int atx,int aty,int size);
void read_initial_transformations(int atx,int aty,int size);
void read_transformations_nonlinear(int atx,int aty,int size);

void writeimage_pgm(char *, PIXEL **, int,int);
void writeimage_raw(char *, PIXEL **, int,int);
void writeimage_pipe(FILE *, PIXEL **, int,int);
void contraction(double **, PIXEL **, int,int);
void smooth_image();
void zooming(double);
void help_dec();
void iterative_decoding(int,int,double);
void iterative_decoding_new_init_image(int,int,double);
void iterative_decoding_testing(int,int,double);
void iterative_decoding_LumInv(int,int,double);
void iterative_decoding_LumInv2(int,int,double);
void iterative_decoding_nonlinear(int,int,double);
void piramidal_decoding(int);
void piramidal_decoding_testing(int);
void piramidal_decoding_LumInv(int);
void piramidal_decoding_LumInv2(int);
void piramidal_decoding_nonlinear(int level);
int quan(double val);
int max_2(int s1, int s2);
int min_2(int s1, int s2);


void read_transformations_adaptive(int ,int ,int,int );
void adaptiveFisherIndexing_2(int ,int);
void adaptiveNewclass(int , int , double **, 
                       int *, int *);
void adaptiveNewclass_2(int , int , double **, 
                       int *, int *);
int adaptiveVariance_class(int x_size,int y_size,double **block);
double variance_3(int x_size,int y_size, double **block, int atx, int aty);
void adaptiveFlips(int x_size,int y_size,double **block,double **flip_block,int iso);

double adaptiveFisherCoding(int,int,int,int, int *,int *,int *, int *,int *);
void traverseImage(int , int , int , int);
void traverseImage_2(int , int , int , int);
void compressRange(int, int , int , int );

void decompressRange(int, int ,int , int);

int getL1Class(int , int, int );

double modified_FisherCoding(int,int,int ,int*,int*,int*,int*,int*);

/* HV patitioning */
void HV_FisherIndexing(int ,int);
double HV_FisherCoding(int ,int ,int ,int ,int* ,int* ,int* ,int* ,int*);
void HV_compressRange(int, int , int , int );
void HV_traverseImage(int , int , int , int);
void HV_traverseImage_2(int , int , int , int);
void HV_decompressRange(int ,int ,int ,int );


void push(struct c** head_ref, int new_data);
void printList(struct node *node);
struct c *getTail(struct c *cur);
struct c *partition(struct c *head, struct c *end,
                       struct c **newHead, struct c **newEnd);
struct c *quickSortRecur(struct c *head, struct c *end);
void quickSort(struct c **headRef);
void FisherIndexing_domainSort(int ,int);
double modified_FisherCoding_1(int atx,int aty,int size,int *xd,int *yd,int *is,
                                                       int* qalf,int *qbet);
double modified_FisherCoding_2(int atx,int aty,int size,int *xd,int *yd,int *is,
                                                       int* qalf,int *qbet);


/* progressive fractal coding */

