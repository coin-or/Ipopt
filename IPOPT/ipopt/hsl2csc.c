/*
  Copyright (C) 2003, International Business Machines and others.
  All Rights Reserved.
  This code is published under the Common Public License.
*/

/* $Id: hsl2csc.c 531 2004-03-11 01:31:07Z andreasw $ */

/* This little function converts the Harwell sparse matrix format into
   symmetric CSC format */

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

void F77_FUNC(hsl2csc,HSL2CSC)
           (fint *n, fint *nnzhsl, fint *irn, fint *jcn,
	    fint *nnzcsc, fint* ia, fint* ja, fint *phsl2csc, fint *ierr)
{
  int i, j, nnzcsc_max, ihsl;
  hslentry *A;

  A = (hslentry*)malloc((*nnzhsl)*sizeof(hslentry));
  if( !A )
    {
      *ierr = -1;
      return;
    }

/* copy data to sort array */
  for (i=0; i<*nnzhsl; i++) {
    if( irn[i]>=jcn[i] )
      {
	A[i].i = irn[i];
	A[i].j = jcn[i];
      }
    else
      {
	A[i].i = jcn[i];
	A[i].j = irn[i];
      }
    A[i].pos = i+1;
  }

  /*  printf("Before qsort:\n\n");
      for(i=0;i<*nnzhsl; i++)
      {
      printf("i = %d  j = %d pos = %d\n",A[i].i,A[i].j,A[i].pos);
      }
  */

  /* sort the HSL entries */
  qsort(A, *nnzhsl, sizeof(hslentry), compare);

  /*  printf("After qsort:\n\n");
      for(i=0;i<*nnzhsl; i++)
      {
      printf("i = %d  j = %d pos = %d\n",A[i].i,A[i].j,A[i].pos);
      }
  */

  nnzcsc_max = *nnzcsc;

  ihsl = 0;
  *nnzcsc = 0;
  for( j=1; j<=*n; j++ )
    {
      ia[j-1] = *nnzcsc + 1;
      /* need a diagonal entry in any case */
      if( *nnzcsc == nnzcsc_max )
      {
	printf("nnzcsc = %d nnzcsc_max = %d nnzhsl = %d\n",*nnzcsc,nnzcsc_max,*nnzhsl);
	/* Not enough space in ia and ja */
	*ierr = -2;
	return;
      }
      ja[*nnzcsc] = j;
      (*nnzcsc)++;

      while( ihsl < *nnzhsl && A[ihsl].j == j)
      {
	if( A[ihsl].i != ja[*nnzcsc-1] )
	  {
	    /* This is a new element */
	    if( *nnzcsc == nnzcsc_max )
	      {
		/* Not enough space in ia and ja */
		printf("nnzcsc = %d nnzcsc_max = %d\n",*nnzcsc,nnzcsc_max);
		*ierr = -2;
		return;
	      }
	    ja[*nnzcsc] = A[ihsl].i;
	    (*nnzcsc)++;
	  }
	phsl2csc[A[ihsl].pos-1] = *nnzcsc;
	ihsl++;
      }
    }
  ia[*n] = *nnzcsc + 1;

  /* That should be it */
  free(A);
  *ierr = 0;
}

#ifdef __cplusplus
}
#endif
