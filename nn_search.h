/************************************************************************
    frap/defrap - fractal image compressor/decompressor
    Copyright (C) 1997 Matthias Ruhl

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
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA

    The author can be contacted at:
      Universitaet Freiburg
      Institut fuer Informatik
      Am Flughafen 17, Geb. 51
      D-79110 Freiburg
      GERMANY
    or via e-mail: ruhl@informatik.uni-freiburg.de
*************************************************************************/

#ifndef NN_SEARCH_H
#define NN_SEARCH_H

#define BUCKETSIZE 10

typedef struct kdtree
{
  int cutdim;
  float cutval;
  float *max,*min;
  struct kdtree *left,*right;
  int *idx,num;
} kdtree;

kdtree *kdtree_build(float **p, int num, int dim);
void kdtree_free(kdtree *t);
int kdtree_search(float *q, float **p, int dim, kdtree *t,
		   float eps, int num, int *nlist);

#endif

