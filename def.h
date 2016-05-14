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



#define IDENTITY      0      /*                          */
#define L_ROTATE90    1      /*                          */
#define R_ROTATE90    2      /*                          */
#define ROTATE180     3      /*       Isometries         */
#define R_VERTICAL    4      /*                          */
#define R_HORIZONTAL  5      /*                          */
#define F_DIAGONAL    6      /*                          */
#define S_DIAGONAL    7      /*                          */

#define matrix_allocate(matrix, hsize, vsize, TYPE) {\
    TYPE *imptr; \
    int _i; \
    matrix = (TYPE **)malloc((vsize)*sizeof(TYPE *));\
    imptr = (TYPE*)malloc((long)(hsize)*(long)(vsize)*sizeof(TYPE));\
    if (imptr == NULL) \
	fatal("\nNo memory in matrix allocate."); \
    for (_i = 0; _i<vsize; ++_i, imptr += hsize) \
	matrix[_i] = imptr; \
 }


#define swap(a,b,TYPE)      {TYPE _temp; _temp=b; b=a; a= _temp;}

#define bound(a)   ((a) < 0.0 ? 0 : ((a)>255.0? 255 : a))

#define MAX_NEIGHBOURS 1000

#define MassCenter     0
#define SaupeFisher    1
#define Saupe          2
#define Fisher         3
#define Hurtgen        4
#define McSaupe        5 

#define TWOPI 6.2831853 

#define VERSION "Mars -- Version 1.0"
#define AUTHOR  "Copyright (C) 1998 Mario Polvere <marpol@iname.com>"



