/*
  Copyright (C) 2003, International Business Machines and others.
  All Rights Reserved.
  This code is published under the Common Public License.
*/

/* $Id: hsl2gcsc.c 531 2004-03-11 01:31:07Z andreasw $ */

/* This little function converts the Harwell sparse matrix format into
   nonsymmetric CSC format */

#ifdef __cplusplus
extern "C" {
#endif

#include <config.h>

#include <stdio.h>
#if STDC_HEADERS
# include <stdlib.h>
# include <stddef.h>
#else
# if HAVE_STDLIB_H
#  include <stdlib.h>
# endif
#endif
#if HAVE_STRING_H
# if !STDC_HEADERS && HAVE_MEMORY_H
#  include <memory.h>
# endif
# include <string.h>
#endif
#if HAVE_STRINGS_H
# include <strings.h>
#endif
#if HAVE_INTTYPES_H
# include <inttypes.h>
#else
# if HAVE_STDINT_H
#  include <stdint.h>
# endif
#endif
#if HAVE_UNISTD_H
# include <unistd.h>
#endif

#include <Ipopt.h>

typedef struct {
    fint i, j, pos;
} hslentry;

static int compare(const void *A, const void *B) {
    hslentry *a, *b;
    a = (hslentry *) A;
    b = (hslentry *) B;
    if (a->j < b->j) return(-1);
    if (a->j > b->j) return(1);
    if (a->i < b->i) return(-1);
    if (a->i > b->i) return(1);
    return(0);
}

void F77_FUNC(hsl2gcsc,HSL2GCSC)
              (fint *n, fint *nnzhsl, fint *irn, fint *jcn,
	       fint *nnzgcsc, fint* ia, fint* ja, fint *phsl2gcsc,
	       fint *phsl2gcsc2, fint *ierr)
{
  int i, j, nnzgcsc_max, nnzgcsc_large, ihsl;
  hslentry *A;

  nnzgcsc_max = *nnzgcsc;
  A = (hslentry*)malloc(2*(*nnzhsl)*sizeof(hslentry));
  if( !A )
    {
      *ierr = -1;
      return;
    }

  ihsl = 0;
/* copy data to sort array */
  for (i=0; i<*nnzhsl; i++) {
    if( irn[i]==jcn[i] )
      {
	A[ihsl].i = irn[i];
	A[ihsl].j = jcn[i];
	A[ihsl].pos = i+1;
	ihsl++;
      }
    else
      {
	A[ihsl].i = jcn[i];
	A[ihsl].j = irn[i];
	A[ihsl].pos = i+1;
	ihsl++;
	A[ihsl].i = irn[i];
	A[ihsl].j = jcn[i];
	A[ihsl].pos = i+1;
	ihsl++;
      }
  }
  nnzgcsc_large = ihsl;

  /*#define DEBUG*/
#ifdef DEBUG
  printf("Before qsort:\n\n");
  for(i=0;i<nnzgcsc_large; i++)
    {
      printf("i = %d  j = %d pos = %d\n",A[i].i,A[i].j,A[i].pos);
    }
#endif

  /* sort the HSL entries */
  qsort(A, nnzgcsc_large, sizeof(hslentry), compare);

#ifdef DEBUG
  printf("After qsort:\n\n");
  for(i=0;i<nnzgcsc_large; i++)
    {
      printf("i = %d  j = %d pos = %d\n",A[i].i,A[i].j,A[i].pos);
    }
#endif

  for( j=0; j<*nnzhsl; j++ )
    {
      phsl2gcsc [j] = -1;
      phsl2gcsc2[j] = -1;
    }

  ihsl = 0;
  *nnzgcsc = 0;
  for( j=1; j<=*n; j++ )
    {
      ia[j-1] = *nnzgcsc + 1;

      while( ihsl < nnzgcsc_large && A[ihsl].j == j)
      {
	if( ihsl == 0 || A[ihsl].i != ja[*nnzgcsc-1] )
	  {
	    /* This is a new element */
	    if( *nnzgcsc == nnzgcsc_max )
	      {
		/* Not enough space in ia and ja */
		printf("nnzgcsc = %d nnzgcsc_max = %d\n",*nnzgcsc,nnzgcsc_max);
		*ierr = -2;
		return;
	      }
	    ja[*nnzgcsc] = A[ihsl].i;
	    (*nnzgcsc)++;
	  }
	i = A[ihsl].pos-1;
	if( phsl2gcsc[i]==-1 )
	  phsl2gcsc [i] = *nnzgcsc;
	else
	  phsl2gcsc2[i] = *nnzgcsc;
	ihsl++;
      }
    }
  ia[*n] = *nnzgcsc + 1;

  /* That should be it */
  free(A);
  *ierr = 0;
}

#ifdef __cplusplus
}
#endif
