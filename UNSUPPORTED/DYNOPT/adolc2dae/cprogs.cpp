/*
  Copyright (C) 2002, Carnegie Mellon University and others.
  All Rights Reserved.
  This code is published under the Common Public License.
*/

/*  $Id: cprogs.cpp 531 2004-03-11 01:31:07Z andreasw $  */
#include "adolc.h"
#include "f2c_adolc.h"

#ifdef AIX
#define myjac_ myjac
#define myjacpat_ myjacpat
#endif

  /*  My Jacobain routine */

extern "C" {

#define maxinc(a,b) if ((a) < (b)) (a) = (b)

integer myjac_(integer* tag,
	       integer* depen,
	       integer* indep,
	       double *argument,
	       integer *jacdim,
	       double *jac)
{
  int TAG = *tag;
  int DEPEN = *depen;
  int INDEP = *indep;
  int JACDIM = *jacdim;

  integer rc= -1;
  double** Jac = myalloc2(DEPEN,INDEP);
  int i, j;

  if( JACDIM < DEPEN )
    {
      fprintf(DIAG_OUT,"Error in myjac: jacdim = %i too small.\n  It must be at least as large as depen = %i.\n", JACDIM, DEPEN);
      exit (-1);
    }
  rc= (integer)jacobian(TAG,DEPEN,INDEP,argument,Jac);
  /* Now copy values into Fortran array */
  for (j=1; j<INDEP; j++)
    {
    for (i=0; i<DEPEN; i++)
      jac[i] = Jac[i][j];
    jac += JACDIM;
    }

  myfree2(Jac);
  return rc;
}

integer myjacpat_(integer* tag,    /* I: tag of tape */
		  integer* depen,  /* I: number of dep vars */
		  integer* indep,  /* I: number of indep vars */
		  integer *nz,     /* I: size of avar and acon; 
				      O: number of nonzero elements */
		  integer *avar,   /* I: array of size nz;  O: sparsity structure */
		  integer *acon )  /* I: array of size nz;  O: sparsity structure */
  /* returns values of ADOLC's jac_pat, or is set to -4 if nz too small */
  /* i.e. there is a problem if return values is negative */
{
  int TAG = *tag;
  int DEPEN = *depen;
  int INDEP = *indep;
  int NZ = *nz;

  integer rc= -1;

  /*
    Call jac_pat to obtain sparsity structure in their format
  */
  double *x = NULL;  /* Since we want the "safe" structure, no point needs
			to be given */

  unsigned int *rb = NULL;  /* We want to consider all variables as
			       single blocks */
  unsigned int *cb = NULL;  /* We want to consider all variables as
			       single blocks */
  unsigned int **crs = new unsigned int* [ DEPEN ]; 
                            /* In here will be the sparsity information */
  int options[2] = {0, 0};

  rc = (integer) jac_pat( TAG, DEPEN, INDEP, x, rb, cb, crs, options );
  if( rc < 0 )
    {
      /* Something went wrong... */
      *nz = 0;
      return rc;
    }

  int nzmax = NZ;
  int icon = 0;
  int inz = 0;
  while( inz < nzmax && icon < DEPEN )
    {
      unsigned int *pcon = crs[icon];
      unsigned int nvar = pcon[0];
      for( unsigned int i = 1; i <= nvar && inz < nzmax; i++ )
	      {
	        avar[inz] = pcon[i] + 1;
	        acon[inz] = icon + 1;
	        inz++;
      	}
      icon ++;
    }
  if( inz == nzmax )
    {
      /*
	There is not enough space to store all nonzero elements
      */
      rc = -4;
      fprintf(DIAG_OUT,"Error in myjacpat: nz = %i is not large enough.\n",nzmax);
    }
  else
    {
      *nz = inz;
    }

  return rc;
}

#ifndef OLDADOLC
/*
  Routine for computing several Hessian-vector products at onces
  (based on Andrea Walter's example hessmat.C)
*/

integer myhessmat_(integer* tag,    /* I: tag of tape */
		   integer* m,      /* I: number of dep vars */
		   integer* n,      /* I: number of indep vars */
		   integer* nrhs,   /* I: number of right hand sides */
		   double *x,   /* I: point at which Hessians are to be evaluated */
		   double *lam, /* I: Lagrange multipliers */
		   double *rhs, /* I: right hand sides (Fortran ordering) */
		   double *res,  /* O: product of Hessians with right hand sides */
		   integer* ldrs /* I: leading dimension of rhs and res */
	       )
  /* Return value:  output of hov_wk_forward and hosv_reverse */
{
  int rc;

  int TAG = *tag;
  int M = *m;
  int N = *n;
  int NRHS = *nrhs;
  int LDRS = *ldrs;

  const int d = 1;
  const int keep = 2;

  int i, j;

  double* yp = new double[M];           /* passive depends (dummy)           */

  double*** Xppp = myalloc(N,NRHS,d);   /* matrix on right-hand side         */
  double*** Zppp = myalloc(NRHS,N,d+1); /* result of Up x H x XPPP           */
  double*** Yppp = myalloc(M,NRHS,d);   /* results of needed hos_wk_forward  */

  /*
    Copy rhs into Xppp
  */
  for (j=0;j<NRHS;j++)
    {
      Xppp[0][j][0] = 0;  /* dummy time variable set to zero */
      for (i=1; i<N; i++)
	Xppp[i][j][0] = rhs[(i-1)+j*LDRS];
    }


  /* The new drivers. First, hov_wk_forward(..) is called.
     So far, it was impossible to store the results of
     a higher-order-vector (=hov) forward in order to perform
     a corresponding reverse sweep (for no particular reason.
     Now we have hov with keep (=wk) and the results needed on
     the way back are stored in a specific tape */

  rc = hov_wk_forward(TAG,M,N,d,keep,NRHS,x,Xppp,yp,Yppp);

  /* The corresponding reverse sweep
     So far we had only a higher-order-scalar (=hos, scalar because
     only one vector on the left-hand-side) for a scalar forward
     call.
     Now, we use the stored vector information (= hos vector)
     to compute multiple lagra_hess_vec at once */

  maxinc(rc, hosv_reverse(TAG,M,N,keep-1,NRHS,lam,Zppp));

  /*
    Copy product from Zppp into res
  */

  for (i=1; i<N; i++)
    for (j=0;j<NRHS;j++)
      res[(i-1)+j*LDRS] = Zppp[j][i][1];

  /*
    Free allocated memory
  */
  myfree3( Yppp );
  myfree3( Zppp );
  myfree3( Xppp );

  delete[] yp;

  return (integer) rc;
}
#endif

}
