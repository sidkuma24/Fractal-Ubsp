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
double Mc_SaupeCoding(int,int,int,int *,int *,int *,int *,int *);
double MassCenterCoding(int,int,int,int *,int *,int *,int *,int *);
double Saupe_FisherCoding(int,int,int,int *,int *,int *,int *,int *);
double entropy(int, int,int,int);
double variance(int, int,int,int);
double variance_2(int, double **, int, int);
int hurtgen_class(int ,double **);
int pack(int , long , FILE *);
int variance_class(int ,double **);
kdtree *kdtree_build(float **,int,int);
void ComputeAverageFactorMc();
void ComputeFeatVectDimSaupe();
void ComputeMcVectors(double **,double **,int,int,double *);
void ComputeSaupeVectors(double **,int,int,float *);
void FisherIndexing(int ,int);
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
void readimage_raw(char *);
void readimage_pgm(char *,int *,int *);
void help_enc();

long unpack(int size, FILE *fin);
void read_transformations(int atx,int aty,int size);
void writeimage_pgm(char *, PIXEL **, int,int);
void writeimage_raw(char *, PIXEL **, int,int);
void writeimage_pipe(FILE *, PIXEL **, int,int);
void contraction(double **, PIXEL **, int,int);
void smooth_image();
void zooming(double);
void help_dec();
void iterative_decoding(int,int,double);
void piramidal_decoding(int);


