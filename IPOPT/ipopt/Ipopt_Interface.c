/*
  Copyright (C) 2003, Kirk Abbott, International Business Machines, and others.
  All Rights Reserved.
  This code is published under the Common Public License.
*/

/* $Id: Ipopt_Interface.c 657 2004-10-07 16:12:51Z andreasw $ */

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

#include "Ipopt_Interface.h"

/* The following determines the length of the Fortran strings */
#define IPOPT_STRLEN 40

/* The following numbers are used if USE_MALLOC is not defined.
   If USE_MALLOC is undefined, IPOPT will try to allocate as much memory
   as possible, with a maximum of IPOPT_MEMMAXTRIAL bytes. Once the maximal
   avaliable memory is determined, that number times IPOPT_MEMFRAC
   is the amount of memory allocated for IPOPT.
*/

#define IPOPT_MEMMAXTRIAL 0x20000000
#define IPOPT_MEMTRIALTOL 0x100000
#define IPOPT_MEMFRAC 0.8

/* Helpful macros for memory allocation */
#define IPOPT_ALLOC(what,number,type,retval) {what = (type *)malloc(number*sizeof(type));if(!what){if(IPOPT_verbose)printf("Ipopt_Interface: Error allocating memory for ""what"" (%d bytes).\n",number*sizeof(type));return retval;}}

#define IPOPT_REALLOC(what,number,type,retval) {what = (type *)realloc(what,number*sizeof(type));if(!what){if(IPOPT_verbose)printf("Ipopt_Interface: Error reallocating memory for ""what"" (%d bytes).\n",number*sizeof(type));return retval;}}

/* Prototype for IPOPT's Fortran subroutine/function calls */
real F77_FUNC(ipopt,IPOPT)
       (fint *N, real *X, fint *M, fint *NLB, fint *ILB,
        real *BNDS_L, fint *NUB, fint *IUB, real *BNDS_U,
        real *V_L, real *V_U, real *LAM, real *C, fint *LRW,
        real *RW, fint *LIW, fint *IW, fint *ITER,
        fint *IERR, void *EV_F(), void *EV_C(), void *EV_G(),
	void *EV_A(), void *EV_H(), void *EV_HLV(), void *EV_HOV(),
	void *EV_HCV(), real *DAT, fint *IDAT,
        fint *NARGS, real *ARGS, char *CARGS, int lCARGS) ;

void F77_FUNC_(ev_hlv_dummy,EV_HLV_DUMMY)
       (fint *task, fint *n, real *x, fint *m, real *lambda,
	real *vin, real *vout, real *DAT, fint *IDAT);

void F77_FUNC_(ev_hov_dummy,EV_HOV_DUMMY)
       (fint *task, fint *n, real *x, fint *m,
	real *vin, real *vout, real *DAT, fint *IDAT);

void F77_FUNC_(ev_hcv_dummy,EV_HCV_DUMMY)
       (fint *task, fint *n, real *x, fint *m, real *lambda,
	real *vin, real *vout, real *DAT, fint *IDAT);

/*
************************************************************************
*  The following is the data structure that has to be set up before a  *
*  call of IPOPT, and will contain the result after a run.             *
************************************************************************
*/
struct _Ipopt {
/* The first 6 entries are set by calling Ipopt_Create */
fint n;                       /* n: number of optimization variables */
fint m;                       /* m: number of (equality) constraints */
fint nlb;                     /* number lower bounds */
fint nub;                     /* number upper bounds */

/* The next 5 entries have to be set by user (note, that memory will be
allocated with the call of Ipopt_Create, so the user does not have
to call malloc again... */
fint *ilb;                    /* indices of lower bounds */
real *bnds_l;                 /* values of lower bounds */
fint *iub;                    /* index of upper bounds */
real *bnds_u;                 /* upper bounds */

/* These function pointers are for the problem statement evalutation
   functions */
pEval_F Eval_F;                /* Function for evaluating objective function */
pEval_C Eval_C;                /* Function for evaluating constraint
                                  functions */
pEval_G Eval_G;                /* Function for evaluating objective gradient */
pEval_A Eval_A;                /* Function for evaluating constraint
                                  Jacobian */
pEval_H Eval_H;                /* Function for evaluating Lagrangian Hessian */

/*
* Option handling stuff (for INIT_PARAMS)
* Use Ipopt_AddParam to set those
*/
char *cargs;                  /* character args */
int  nargs_alloc;             /* currently allocated length for cargs
                                 and args */
fint nargs;                   /* args length */
real *args;                   /* double args */
};

/**
 * funtion for copying C-string into fortran character array.
 * Taken from ipoptAMPL.c
 */
static void str2fstr(char *c_str, char *f_str, size_t ftnlen)
{
  size_t i, len;

  len = strlen(c_str);
  if (len > ftnlen) {
    len = ftnlen;
  }
  memcpy(f_str, c_str, len);
  for (i = len; i < ftnlen; f_str[i++] = ' ');
}

static int IPOPT_verbose = 1;
void DLLEXPORT Ipopt_SetVerbose(int verbose)
{
  IPOPT_verbose = verbose;
}

/**
 * Ipopt_Version.
 * Returns the version of this code
 */
char DLLEXPORT *Ipopt_Version()
{
  static char *IPOPTversion = PACKAGE_VERSION;
  return IPOPTversion;
}

/**
 * Ipopt_Create.
 * Creates the problem.
 * The philosophy here is that all memory that is needed should
 * be allocated at the end of this call. The values will be filled
 * in later.
 *
 * @param int n, number of variables n
 * @param int m, number of constraints m
 * @param int nlb, number of lower bounds
 * @param int nub, number of upper bounds
 */
Ipopt DLLEXPORT Ipopt_Create(fint n, fint m, fint nlb,
			     fint *ilb, real *bnds_l, fint nub,
			     fint *iub, real *bnds_u, pEval_F Eval_F,
			     pEval_C Eval_C, pEval_G Eval_G,
			     pEval_A Eval_A, pEval_H Eval_H)
{
  Ipopt result;
  int i;

  /*result = (Ipopt)malloc(sizeof(struct Ipopt));*/
  IPOPT_ALLOC(result,1,struct _Ipopt,NULL);

  if( n<1 || m<0 || nlb<0 || nub<0 ) {
    printf("Ipopt_Interface: Invalid input for Ipopt_Create.\n");
    return NULL;
  }
  result->n = n;
  result->m = m;
  result->nlb = nlb;
  result->nub = nub;

  if( nlb > 0 ) {
    IPOPT_ALLOC(result->ilb,nlb,fint,NULL);
    IPOPT_ALLOC(result->bnds_l,nlb,double,NULL);
    for( i=0; i<nlb; i++) {
      result->bnds_l[i] = bnds_l[i];
      result->ilb[i] = ilb[i];
    }
  }

  if( nub > 0 ) {
    IPOPT_ALLOC(result->iub,nub,fint,NULL);
    IPOPT_ALLOC(result->bnds_u,nub,double,NULL);
    for( i=0; i<nub; i++) {
      result->bnds_u[i] = bnds_u[i];
      result->iub[i] = iub[i];
    }
  }

  result->nargs = 0;
  result->nargs_alloc = 0;

  /* Set the function pointers for the evaluation functions */
  result->Eval_F=Eval_F;
  result->Eval_C=Eval_C;
  result->Eval_G=Eval_G;
  result->Eval_A=Eval_A;
  result->Eval_H=Eval_H;

  return result;
}

/**
 * Ipopt_Destroy.
 * Destroys the problem structure.
 * @param Ipopt problem, Ipopt instance to be deleted
 */
void Ipopt_Destroy(Ipopt problem)
{
  if (problem != NULL) {
    if (problem->nlb > 0) {
      free(problem->ilb);
      free(problem->bnds_l);
    }

    if (problem->nub > 0) {
      free(problem->iub);
      free(problem->bnds_u);
    }

    if (problem->nargs_alloc > 0) {
      free(problem->args);
      free(problem->cargs);
    }

    free(problem);
  }
}

/**
 * IPOPT_EVAL_F
 * Computation of objective function
 */
void EV_F(fint *n, real *x, real *f, real *DAT, fint *IDAT)
{
  int status;
  void *userData = (void*)DAT;
  Ipopt problem=(Ipopt)IDAT;
  status = (problem->Eval_F)(*n, x, f, userData);
  if( status ) {
    printf("Warning: Eval_F returned status = %d.\n",status);
  }
}

/**
 * IPOPT_EVAL_G
 * Computation of gradient of objective function
 */
void EV_G(fint *n, real *x, real *g, real *DAT, fint *IDAT)
{
  int status;
  void *userData = (void*)DAT;
  Ipopt problem=(Ipopt)IDAT;
  status = (problem->Eval_G)(*n, x, g, userData);
  if( status ) {
    printf("Warning: Eval_G returned status = %d.\n",status);
  }
}

/**
 * IPOPT_EVAL_C
 * Computation of equality constraints
 */
void EV_C(fint *n, real *x, fint *m, real *c, real *DAT, fint *IDAT)
{
  int status;
  void *userData = (void*)DAT;
  Ipopt problem=(Ipopt)IDAT;
  status = (problem->Eval_C)(*n, x, *m, c, userData);
  if( status ) {
    printf("Warning: Eval_C returned status = %d.\n",status);
  }
}

/**
 * IPOPT_EVAL_A
 * Computation of Jacobian of equality constraints
 */
void EV_A(fint *task, fint *n, real *x, fint *nz,
	  real *A, fint *Acon, fint *Avar, real *DAT, fint *IDAT)
{
  int status;
  void *userData = (void*)DAT;
  Ipopt problem=(Ipopt)IDAT;
  status = (problem->Eval_A)(*task, *n, x, nz, A, Acon, Avar, userData);
  if( status ) {
    printf("Warning: Eval_A returned status = %d.\n",status);
  }
}

/**
 * IPOPT_EVAL_H
 * Computation of Hessian of Lagrangian
 */
void EV_H(fint *task, fint *n, real *x, fint *m,
	  real *lambda, fint *nnzh, real *hess,
	  fint *irnh, fint *icnh, real *DAT, fint *IDAT)
{
  int status;
  void *userData = (void*)DAT;
  Ipopt problem=(Ipopt)IDAT;

  if( problem->Eval_H == NULL ) {
    printf("Error: No Eval_H function given to Ipopt_Create.  In that case, you need to\n       choose options that do not require second derivatives.\n");
    exit(1);
  }
  status = (problem->Eval_H)(*task, *n, x, *m, lambda, nnzh, hess, irnh, icnh,
		             userData);
  if( status ) {
    printf("Warning: Eval_H returned status = %d.\n",status);
  }
}

void EV_HLV(fint *task, fint *n, real *x, fint *m, real *lambda,
	    real *vin, real *vout, real *DAT, fint *IDAT)
{
  F77_FUNC_(ev_hlv_dummy,EV_HLV_DUMMY)
    (task, n, x, m, lambda, vin, vout, DAT, IDAT);
}

void EV_HOV(fint *task, fint *n, real *x, fint *m,
	    real *vin, real *vout, real *DAT, fint *IDAT)
{
  F77_FUNC_(ev_hov_dummy,EV_HOV_DUMMY)
    (task, n, x, m, vin, vout, DAT, IDAT);
}

void EV_HCV(fint *task, fint *n, real *x, fint *m, real *lambda,
	    real *vin, real *vout, real *DAT, fint *IDAT)
{
  F77_FUNC_(ev_hcv_dummy,EV_HCV_DUMMY)
    (task, n, x, m, lambda, vin, vout, DAT, IDAT);
}

/**
 * Ipopt_Solve.
 * Solves the provided problem.
 * @param Ipopt *problem
 * @param void *userData
 * @param  real *Obj, final value of objective function
 * @return int, the solution status
 */
int Ipopt_Solve(Ipopt problem, void *userData, real *x,
		real *Obj, real *v_lb, real *v_ub,
		real *lambda, real *c, fint *iter)
{
  fint *iw, liw, lrw, idummy, retval;
  double *rw;

#ifndef USE_MALLOC
  int memsize, memmin, memmax;
  void *ptr;
  /* Try to allocate almost as much memory as possible */
  memmin = 0x0;
  memmax = IPOPT_MEMMAXTRIAL;
  while( memmax - memmin > IPOPT_MEMTRIALTOL )
    {
      memsize = (memmax + memmin)/2;
      ptr = malloc(memsize);
      if( ptr )
	{
	  memmin = memsize;
	  free( ptr );
	}
      else
	{
	  memmax = memsize;
	}
    }
  memsize = (int)( ((double)memmin)*IPOPT_MEMFRAC );
  liw = memsize/(sizeof(fint)+sizeof(real));
  lrw = liw;

  /* Allocate memory for the work space */
  IPOPT_ALLOC(iw,liw,fint,-1);
  IPOPT_ALLOC(rw,lrw,real,-2);
#else /* USE_MALLOC */
  liw = 0;
  lrw = 0;
  iw = NULL;
  rw = NULL;
#endif /* USE_MALLOC */

  *Obj = F77_FUNC(ipopt,IPOPT)
              (&(problem->n), x, &(problem->m),
	       &(problem->nlb), problem->ilb, problem->bnds_l,
	       &(problem->nub), problem->iub, problem->bnds_u,
	       v_lb, v_ub, lambda, c, &lrw, rw, &liw, iw, iter,
	       &retval, (void *(*)())EV_F, (void *(*)())EV_C,
	       (void *(*)())EV_G, (void *(*)())EV_A, (void *(*)())EV_H,
	       (void *(*)())EV_HLV, (void *(*)())EV_HOV, (void *(*)())EV_HCV,
	       (real*)userData, (fint*)problem,
	       &(problem->nargs), problem->args, problem->cargs,
	       IPOPT_STRLEN);

  return retval;
}

/**
 * Add a parameter to the list in problem
 * Returns 0 if everything went fine, or returns -1 if there was a problem
 * in allocating memory
 */
int DLLEXPORT Ipopt_AddParam(Ipopt problem, char *name, double value)
{

  if( problem->nargs_alloc == 0 )
    {
      problem->nargs_alloc = 10;
      IPOPT_ALLOC(problem->args,problem->nargs_alloc,real,-1);
      IPOPT_ALLOC(problem->cargs,(IPOPT_STRLEN)*(problem->nargs_alloc),char,-1);
    }
  else if( problem->nargs == problem->nargs_alloc )
    {
      problem->nargs_alloc += 10;
      IPOPT_REALLOC(problem->args,problem->nargs_alloc,real,-1);
      IPOPT_REALLOC(problem->cargs,(IPOPT_STRLEN)*(problem->nargs_alloc),char,-1);
    }

  (problem->args)[problem->nargs] = (real)value;
  str2fstr(name,&((problem->cargs)[IPOPT_STRLEN*problem->nargs]),IPOPT_STRLEN);
  problem->nargs++;

  return 0;
}


#ifdef __cplusplus
}
#endif
