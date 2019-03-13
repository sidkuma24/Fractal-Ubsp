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

#include <stdlib.h>
#include <assert.h>
#include "nn_search.h"

#include "def.h"

#define HEAPMAX 10000

typedef struct
{
  float dist;
  kdtree *t;
} heapT;

static heapT *h = NULL;
static int hsize;
static int cutdim;
static float **p;

/* Heap-functions for kdtree_search */

static void hinit()
{
  if(h == NULL)
    h = (heapT *) malloc(sizeof(heapT)*HEAPMAX);
  hsize = 0;
}

static void hfree()
{
  if(h != NULL) free(h);
}

static int hempty()
{
  return (hsize == 0);
}

static void hpush(kdtree *t, float *q, int dim)
{
  int i,pos;
  float d,tmp;

  assert(hsize < HEAPMAX);
  d = 0.0;
  for(i=0;i<dim;i++)
    {
      if(q[i] < t->min[i]) 
	{
	  tmp = t->min[i]-q[i];
	  d += tmp*tmp;
	}
      else
	if(q[i] > t->max[i]) 
	  {
	    tmp = q[i]-t->max[i];
	    d += tmp*tmp;
	  }
    }
  pos = ++hsize;
  while(pos != 1 && h[pos/2].dist > d)
    {
      h[pos] = h[pos/2];
      pos /= 2;
    }
  h[pos].dist = d;
  h[pos].t = t;
}

static kdtree *hpop()
{
  int pos;
  kdtree *tmp;

  assert(hsize != 0);
  tmp = h[1].t;
  if(--hsize != 0) 
    {
      pos = 1;
      while(2*pos <= hsize)
	{
	  if(2*pos+1 <= hsize && h[2*pos+1].dist < h[2*pos].dist)
	    if(h[2*pos+1].dist < h[hsize+1].dist)
	      {
		h[pos] = h[2*pos+1]; pos = 2*pos+1; continue;
	      }
	    else
	      break;
	  if(h[2*pos].dist < h[hsize+1].dist)
	    {
	      h[pos] = h[2*pos]; pos *= 2; continue;
	    }
	  break;
	}
      h[pos] = h[hsize+1];
    }
  return tmp;
}

/* Find 'num' (1+'eps')-nearest neighbors in kdtree 't' (built from
   points 'p', and put them in 'nlist' */

int kdtree_search(float *q, float **p, int dim, kdtree *t,
		   float eps, int num, int *nlist)
{
  int i,j,found;
  float d, tmp,*ptmp,*qtmp,eps2;
  float dist2[MAX_NEIGHBOURS];

  hinit();
  hpush(t,q,dim);
  found = 0;
  eps2 = (1.0+eps)*(1.0+eps);
  while(!hempty() && (found < num || eps2*h[1].dist < dist2[num-1]))
    {
      t = hpop();
      while(t->num == 0)
	{
	  if(q[t->cutdim] < t->cutval)
	    {
	      hpush(t->right,q,dim);
	      t = t->left;
	    }
	  else
	    {
	      hpush(t->left,q,dim);
	      t = t->right;
	    }
	}

      /* Leaf node: process members */
      for(i=0;i<t->num;i++)
	{
	  d = 0.0;
	  for(qtmp=q,ptmp=p[t->idx[i]],j=0;j<dim;j++)
	    {
	      tmp = (*qtmp++)-(*ptmp++);
	      d += tmp*tmp;
	    }
	  if(found == num && d > dist2[num-1]) continue;
	  if(found < num) found++;
	  for(j=found-1;j>0 && dist2[j-1]>d;j--)
	    { 
	      dist2[j] = dist2[j-1]; nlist[j] = nlist[j-1]; 
	    }
	  dist2[j] = d;
	  nlist[j] = t->idx[i];
	}
    }
  return found;
}

/* Functions for building a kd-tree */

static int compare(const void *a, const void *b)
{
  float p1,p2;

  p1 = p[*((int *)a)][cutdim];
  p2 = p[*((int *)b)][cutdim];
  return p1<p2?-1:p1>p2;
}

static kdtree *kdtree_build2(int *idx, int num, int dim)
{
  int i,d;
  float spread,min,max;
  kdtree *tmp;

  if(num == 0) return NULL;

  tmp = (kdtree *) malloc(sizeof(kdtree));

  /* 1. Find dimension with greatest spread */

  tmp->min = (float *) malloc(sizeof(float)*dim);
  tmp->max = (float *) malloc(sizeof(float)*dim);

  spread = -1.0;
  for(d=0;d<dim;d++)
    {
      min = max = p[idx[0]][d];
      for(i=1;i<num;i++)
	{
	  if(min > p[idx[i]][d]) min = p[idx[i]][d];
	  else
	    if(max < p[idx[i]][d]) max = p[idx[i]][d];
	}
      if(max-min > spread)
	{
	  spread = max-min;
	  cutdim = d;
	}
      tmp->min[d] = min;
      tmp->max[d] = max;
    }
  
  /* 2a. If the number of points is less then BUCKETSIZE:
        pack all points in the bucket, and we're done */

  if(num <= BUCKETSIZE)
    {
      tmp->cutdim = -1;
      tmp->left = tmp->right = NULL;
      tmp->idx = (int *) malloc(sizeof(int)*num);
      tmp->num = num;
      for(i=0;i<num;i++) tmp->idx[i] = idx[i];
      return tmp;
    }

  /* 2b. Otherwise, we have to build a tree recursively */

  /* 3. Sort vertices according to coordinate 'cutdim' */

  qsort(idx,num,sizeof(int),compare);
  
  /* 4. Split vertices evenly and recursively */

  tmp->cutdim = cutdim;
  tmp->cutval = (p[idx[num/2-1]][cutdim] + p[idx[num/2]][cutdim])/2.0;
  tmp->idx = NULL;
  tmp->num = 0;

  tmp->left  = kdtree_build2(idx,num/2,dim);
  tmp->right = kdtree_build2(idx+num/2,num-num/2,dim);

  return tmp;
}

kdtree *kdtree_build(float **p2, int num, int dim)
{
  int *idx;
  int i;
  kdtree *tmp;

  p = p2;
  idx = (int *) malloc(sizeof(int)*num);
  for(i=0;i<num;i++) idx[i] = i;
  tmp = kdtree_build2(idx,num,dim);
  free(idx);
  return tmp;
}

/* To free the memory allocated by a kd-tree */

void kdtree_free(kdtree *t)
{
  if(t == NULL) return;
  kdtree_free(t->left);
  kdtree_free(t->right);
  if(t->idx) free(t->idx);
  free(t->max);
  free(t->min);
  free(t);
  hfree();
}


