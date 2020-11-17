/* Copyright (C) 2004, 2010 International Business Machines and others.
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * Authors:  Carl Laird, Andreas Waechter     IBM    2004-09-03
 */

#include "IpStdCInterface.h"
#include "IpoptConfig.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* if we build without Fortran compiler around, then make up a Fortran name mangling scheme that often works */
#ifndef F77_FUNC
#define F77_FUNC(name,NAME) name ## _
#endif

/* ToDo: The following needs to be adapted based on configuration */
typedef IPOPT_FORTRAN_INTEGER_TYPE fint;
typedef double fdouble;
typedef float ffloat;
/* in configure, we checked whether an int* is 32 or 64bit long to decide how much space
 * is needed to store a pointer
 * thus, we can use a void* here to represent a pointer for Fortran
 */
typedef void* fptr;

/** Return value for indicating that evaluation could be done without problem. */
static const fint OKRetVal = 0;
/** Return value for indicating that evaluation could not be done without problem. */
static const fint NotOKRetVal = 1;

/* Function pointer types for the Fortran callback functions */
#ifdef IPOPT_SINGLE
typedef void (*FEval_F_CB)(
   fint*    N,
   ffloat*  X,
   fint*    NEW_X,
   ffloat*  OBJVAL,
   fint*    IDAT,
   ffloat*  DDAT,
   fint*    IERR
);

typedef void (*FEval_G_CB)(
   fint*    N,
   ffloat*  X,
   fint*    NEW_X,
   fint*    M,
   ffloat*  G,
   fint*    IDAT,
   ffloat*  DDAT,
   fint*    IERR
);

typedef void (*FEval_Grad_F_CB)(
   fint*    N,
   ffloat*  X,
   fint*    NEW_X,
   ffloat*  GRAD,
   fint*    IDAT,
   ffloat*  DDAT,
   fint*    IERR
);

typedef void (*FEval_Jac_G_CB)(
   fint*    TASK,
   fint*    N,
   ffloat*  X,
   fint*    NEW_X,
   fint*    M,
   fint*    NNZJAC,
   fint*    IROW,
   fint*    JCOL,
   ffloat*  VALUES,
   fint*    IDAT,
   ffloat*  DDAT,
   fint*    IERR
);

typedef void (*FEval_Hess_CB)(
   fint*    TASK,
   fint*    N,
   ffloat*  X,
   fint*    NEW_X,
   ffloat*  OBJFACT,
   fint*    M,
   ffloat*  LAMBDA,
   fint*    NEW_LAM,
   fint*    NNZHESS,
   fint*    IROW,
   fint*    JCOL,
   ffloat*  VALUES,
   fint*    IDAT,
   ffloat*  DDAT,
   fint*    IERR
);

typedef void (*FIntermediate_CB)(
   fint*    ALG_MODE,
   fint*    ITER_COUNT,
   ffloat*  OBJVAL,
   ffloat*  INF_PR,
   ffloat*  INF_DU,
   ffloat*  MU,
   ffloat*  DNORM,
   ffloat*  REGU_SIZE,
   ffloat*  ALPHA_DU,
   ffloat*  ALPHA_PR,
   fint*    LS_TRIAL,
   fint*    IDAT,
   ffloat*  DDAT,
   fint*    ISTOP
);

struct _FUserData
{
   fint*            IDAT;
   ffloat*          DDAT;
   FEval_F_CB       EVAL_F;
   FEval_G_CB       EVAL_G;
   FEval_Grad_F_CB  EVAL_GRAD_F;
   FEval_Jac_G_CB   EVAL_JAC_G;
   FEval_Hess_CB    EVAL_HESS;
   FIntermediate_CB INTERMEDIATE_CB;
   IpoptProblem     Problem;
};
#else
typedef void (*FEval_F_CB)(
   fint*    N,
   fdouble* X,
   fint*    NEW_X,
   fdouble* OBJVAL,
   fint*    IDAT,
   fdouble* DDAT,
   fint*    IERR
);

typedef void (*FEval_G_CB)(
   fint*    N,
   fdouble* X,
   fint*    NEW_X,
   fint*    M,
   fdouble* G,
   fint*    IDAT,
   fdouble* DDAT,
   fint*    IERR
);

typedef void (*FEval_Grad_F_CB)(
   fint*    N,
   fdouble* X,
   fint*    NEW_X,
   fdouble* GRAD,
   fint*    IDAT,
   fdouble* DDAT,
   fint*    IERR
);

typedef void (*FEval_Jac_G_CB)(
   fint*    TASK,
   fint*    N,
   fdouble* X,
   fint*    NEW_X,
   fint*    M,
   fint*    NNZJAC,
   fint*    IROW,
   fint*    JCOL,
   fdouble* VALUES,
   fint*    IDAT,
   fdouble* DDAT,
   fint*    IERR
);

typedef void (*FEval_Hess_CB)(
   fint*    TASK,
   fint*    N,
   fdouble* X,
   fint*    NEW_X,
   fdouble* OBJFACT,
   fint*    M,
   fdouble* LAMBDA,
   fint*    NEW_LAM,
   fint*    NNZHESS,
   fint*    IROW,
   fint*    JCOL,
   fdouble* VALUES,
   fint*    IDAT,
   fdouble* DDAT,
   fint*    IERR
);

typedef void (*FIntermediate_CB)(
   fint*    ALG_MODE,
   fint*    ITER_COUNT,
   fdouble* OBJVAL,
   fdouble* INF_PR,
   fdouble* INF_DU,
   fdouble* MU,
   fdouble* DNORM,
   fdouble* REGU_SIZE,
   fdouble* ALPHA_DU,
   fdouble* ALPHA_PR,
   fint*    LS_TRIAL,
   fint*    IDAT,
   fdouble* DDAT,
   fint*    ISTOP
);

struct _FUserData
{
   fint*            IDAT;
   fdouble*         DDAT;
   FEval_F_CB       EVAL_F;
   FEval_G_CB       EVAL_G;
   FEval_Grad_F_CB  EVAL_GRAD_F;
   FEval_Jac_G_CB   EVAL_JAC_G;
   FEval_Hess_CB    EVAL_HESS;
   FIntermediate_CB INTERMEDIATE_CB;
   IpoptProblem     Problem;
};
#endif
typedef struct _FUserData FUserData;

static Bool eval_f(
   Index       n,
   Number*     x,
   Bool        new_x,
   Number*     obj_value,
   UserDataPtr user_data
)
{
   fint N = n;
   fint NEW_X = new_x;
   FUserData* fuser_data = (FUserData*) user_data;
   fint* IDAT = fuser_data->IDAT;
#ifdef IPOPT_SINGLE
   ffloat*  DDAT = fuser_data->DDAT;
#else
   fdouble* DDAT = fuser_data->DDAT;
#endif
   fint IERR = 0;

   fuser_data->EVAL_F(&N, x, &NEW_X, obj_value, IDAT, DDAT, &IERR);

   return (Bool) (IERR == OKRetVal);
}

static Bool eval_grad_f(
   Index       n,
   Number*     x,
   Bool        new_x,
   Number*     grad_f,
   UserDataPtr user_data
)
{
   fint N = n;
   fint NEW_X = new_x;
   FUserData* fuser_data = (FUserData*) user_data;
   fint* IDAT = fuser_data->IDAT;
#ifdef IPOPT_SINGLE
   ffloat*  DDAT = fuser_data->DDAT;
#else
   fdouble* DDAT = fuser_data->DDAT;
#endif
   fint IERR = 0;

   fuser_data->EVAL_GRAD_F(&N, x, &NEW_X, grad_f, IDAT, DDAT, &IERR);

   return (Bool) (IERR == OKRetVal);
}

static Bool eval_g(
   Index       n,
   Number*     x,
   Bool        new_x,
   Index       m,
   Number*     g,
   UserDataPtr user_data
)
{
   fint N = n;
   fint NEW_X = new_x;
   fint M = m;
   FUserData* fuser_data = (FUserData*) user_data;
   fint* IDAT = fuser_data->IDAT;
#ifdef IPOPT_SINGLE
   ffloat*  DDAT = fuser_data->DDAT;
#else
   fdouble* DDAT = fuser_data->DDAT;
#endif
   fint IERR = 0;

   fuser_data->EVAL_G(&N, x, &NEW_X, &M, g, IDAT, DDAT, &IERR);

   return (Bool) (IERR == OKRetVal);
}

static Bool eval_jac_g(
   Index       n,
   Number*     x,
   Bool        new_x,
   Index       m,
   Index       nele_jac,
   Index*      iRow,
   Index*      jCol,
   Number*     values,
   UserDataPtr user_data
)
{
   fint N = n;
   fint NEW_X = new_x;
   fint M = m;
   fint NNZJAC = nele_jac;
   fint TASK;
   FUserData* fuser_data = (FUserData*) user_data;
   fint* IDAT = fuser_data->IDAT;
#ifdef IPOPT_SINGLE
   ffloat*  DDAT = fuser_data->DDAT;
#else
   fdouble* DDAT = fuser_data->DDAT;
#endif
   fint IERR = 0;

   if( iRow != NULL && jCol != NULL && values == NULL )
   {
      /* Only request the structure */
      TASK = 0;
   }
   else if( iRow == NULL && jCol == NULL && values != NULL )
   {
      /* Only request the values */
      TASK = 1;
   }
   else
   {
      printf("Error in IpStdFInterface eval_jac_g!\n");
      return (Bool) 0;
   }

   fuser_data->EVAL_JAC_G(&TASK, &N, x, &NEW_X, &M, &NNZJAC, iRow, jCol, values, IDAT, DDAT, &IERR);

   return (Bool) (IERR == OKRetVal);
}

static Bool eval_h(
   Index       n,
   Number*     x,
   Bool        new_x,
   Number      obj_factor,
   Index       m,
   Number*     lambda,
   Bool        new_lambda,
   Index       nele_hess,
   Index*      iRow,
   Index*      jCol,
   Number*     values,
   UserDataPtr user_data
)
{
   fint N = n;
   fint NEW_X = new_x;
   fint M = m;
   fint NEW_LAM = new_lambda;
   fint NNZHESS = nele_hess;
   fint TASK;
   FUserData* fuser_data = (FUserData*) user_data;
   fint* IDAT = fuser_data->IDAT;
#ifdef IPOPT_SINGLE
   ffloat*  DDAT = fuser_data->DDAT;
#else
   fdouble* DDAT = fuser_data->DDAT;
#endif
   fint IERR = 0;

   if( iRow != NULL && jCol != NULL && values == NULL )
   {
      /* Only request the structure */
      TASK = 0;
   }
   else if( iRow == NULL && jCol == NULL && values != NULL )
   {
      /* Only request the values */
      TASK = 1;
   }
   else
   {
      printf("Error in IpStdFInterface eval_hess!\n");
      return (Bool) 0;
   }

   fuser_data->EVAL_HESS(&TASK, &N, x, &NEW_X, &obj_factor, &M, lambda, &NEW_LAM, &NNZHESS, iRow, jCol, values, IDAT, DDAT, &IERR);

   return (Bool) (IERR == OKRetVal);
}

static Bool intermediate_cb(
   Index       alg_mod,
   Index       iter_count,
   Number      obj_value,
   Number      inf_pr,
   Number      inf_du,
   Number      mu,
   Number      d_norm,
   Number      regularization_size,
   Number      alpha_du,
   Number      alpha_pr,
   Index       ls_trials,
   UserDataPtr user_data
)
{
   FUserData* fuser_data = (FUserData*) user_data;
   fint ALG_MODE = alg_mod;
   fint ITER_COUNT = iter_count;
#ifdef IPOPT_SINGLE
   ffloat  OBJVAL = obj_value;
   ffloat  INF_PR = inf_pr;
   ffloat  INF_DU = inf_du;
   ffloat  MU = mu;
   ffloat  DNORM = d_norm;
   ffloat  REGU_SIZE = regularization_size;
   ffloat  ALPHA_DU = alpha_du;
   ffloat  ALPHA_PR = alpha_pr;
#else
   fdouble OBJVAL = obj_value;
   fdouble INF_PR = inf_pr;
   fdouble INF_DU = inf_du;
   fdouble MU = mu;
   fdouble DNORM = d_norm;
   fdouble REGU_SIZE = regularization_size;
   fdouble ALPHA_DU = alpha_du;
   fdouble ALPHA_PR = alpha_pr;
#endif
   fint LS_TRIAL = ls_trials;
   fint* IDAT = fuser_data->IDAT;
#ifdef IPOPT_SINGLE
   ffloat*  DDAT = fuser_data->DDAT;
#else
   fdouble* DDAT = fuser_data->DDAT;
#endif
   fint ISTOP = 0;

   if( !fuser_data->INTERMEDIATE_CB )
   {
      return (Bool) TRUE;
   }

   fuser_data->INTERMEDIATE_CB(&ALG_MODE, &ITER_COUNT, &OBJVAL, &INF_PR, &INF_DU, &MU, &DNORM, &REGU_SIZE, &ALPHA_DU,
                               &ALPHA_PR, &LS_TRIAL, IDAT, DDAT, &ISTOP);

   return (Bool) (ISTOP == OKRetVal);
}

IPOPTLIB_EXPORT fptr F77_FUNC(ipcreate, IPCREATE)(
   fint*           N,
#ifdef IPOPT_SINGLE
   ffloat*         X_L,
   ffloat*         X_U,
#else
   fdouble*        X_L,
   fdouble*        X_U,
#endif
   fint*           M,
#ifdef IPOPT_SINGLE
   ffloat*         G_L,
   ffloat*         G_U,
#else
   fdouble*        G_L,
   fdouble*        G_U,
#endif
   fint*           NELE_JAC,
   fint*           NELE_HESS,
   fint*           IDX_STY,
   FEval_F_CB      EVAL_F,
   FEval_G_CB      EVAL_G,
   FEval_Grad_F_CB EVAL_GRAD_F,
   FEval_Jac_G_CB  EVAL_JAC_G,
   FEval_Hess_CB   EVAL_HESS
)
{
   Index n = *N;
   Index m = *M;
   Index nele_jac = *NELE_JAC;
   Index nele_hess = *NELE_HESS;
   Index index_style = *IDX_STY;

   FUserData* fuser_data;

   fuser_data = (FUserData*) malloc(sizeof(FUserData));

   /* First create a new IpoptProblem object; if that fails return 0 */
   fuser_data->Problem = CreateIpoptProblem(n, X_L, X_U, m, G_L, G_U, nele_jac, nele_hess,
                         index_style, eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h);
   if( fuser_data->Problem == NULL )
   {
      free(fuser_data);
      return (fptr)NULL;
   }

   /* Store the information for the callback function */
   fuser_data->EVAL_F = EVAL_F;
   fuser_data->EVAL_G = EVAL_G;
   fuser_data->EVAL_GRAD_F = EVAL_GRAD_F;
   fuser_data->EVAL_JAC_G = EVAL_JAC_G;
   fuser_data->EVAL_HESS = EVAL_HESS;
   fuser_data->INTERMEDIATE_CB = NULL;

   return (fptr)fuser_data;
}

IPOPTLIB_EXPORT void F77_FUNC(ipfree, IPFREE)(
   fptr* FProblem
)
{
   FUserData* fuser_data = (FUserData*) *FProblem;

   FreeIpoptProblem(fuser_data->Problem);
   free(fuser_data);

   *FProblem = (fptr)NULL;
}

IPOPTLIB_EXPORT fint F77_FUNC(ipsolve, IPSOLVE)(
   fptr*    FProblem,
#ifdef IPOPT_SINGLE
   ffloat*  X,
   ffloat*  G,
   ffloat*  OBJ_VAL,
   ffloat*  MULT_G,
   ffloat*  MULT_X_L,
   ffloat*  MULT_X_U,
#else
   fdouble* X,
   fdouble* G,
   fdouble* OBJ_VAL,
   fdouble* MULT_G,
   fdouble* MULT_X_L,
   fdouble* MULT_X_U,
#endif
   fint*    IDAT,
#ifdef IPOPT_SINGLE
   ffloat*  DDAT
#else
   fdouble* DDAT
#endif
)
{
   FUserData* fuser_data = (FUserData*) *FProblem;
   UserDataPtr user_data;

   fuser_data->IDAT = IDAT;
   fuser_data->DDAT = DDAT;
   user_data = (UserDataPtr) fuser_data;

   return (fint)IpoptSolve(fuser_data->Problem, X, G, OBJ_VAL, MULT_G, MULT_X_L, MULT_X_U, user_data);
}

static char* f2cstr(
   char* FSTR,
   int   slen
)
{
   int len;
   char* cstr;
   for( len = slen; len > 0; --len )
   {
      if( FSTR[len - 1] != ' ' )
      {
         break;
      }
   }
   cstr = (char*) malloc(sizeof(char) * (len + 1));
   strncpy(cstr, FSTR, len);
   cstr[len] = '\0';

   return cstr;
}

/* ToDo make sure position of vlen and klen are at the right place */
IPOPTLIB_EXPORT fint F77_FUNC(ipaddstroption, IPADDSTROPTION)(
   fptr* FProblem,
   char* KEYWORD,
   char* VALUE,
   int   klen,
   int   vlen
)
{
   char* keyword;
   char* val;
   FUserData* fuser_data = (FUserData*) *FProblem;
   fint retval;

   keyword = f2cstr(KEYWORD, klen);
   val = f2cstr(VALUE, vlen);

   retval = AddIpoptStrOption(fuser_data->Problem, keyword, val);

   free(val);
   free(keyword);

   return retval ? OKRetVal : NotOKRetVal;
}

IPOPTLIB_EXPORT fint F77_FUNC(ipaddnumoption, IPADDNUMOPTION)(
   fptr*    FProblem,
   char*    KEYWORD,
#ifdef IPOPT_SINGLE
   ffloat*  VALUE,
#else
   fdouble* VALUE,
#endif
   int      klen
)
{
   char* keyword;
   FUserData* fuser_data = (FUserData*) *FProblem;
   fint retval;

   keyword = f2cstr(KEYWORD, klen);

   retval = AddIpoptNumOption(fuser_data->Problem, keyword, *VALUE);

   free(keyword);

   return retval ? OKRetVal : NotOKRetVal;
}

IPOPTLIB_EXPORT fint F77_FUNC(ipaddintoption, IPADDINTOPTION)(
   fptr* FProblem,
   char* KEYWORD,
   fint* VALUE,
   int   klen
)
{
   char* keyword;
   FUserData* fuser_data = (FUserData*) *FProblem;
   Int value = *VALUE;
   fint retval;

   keyword = f2cstr(KEYWORD, klen);

   retval = AddIpoptIntOption(fuser_data->Problem, keyword, value);

   free(keyword);

   return retval ? OKRetVal : NotOKRetVal;
}

IPOPTLIB_EXPORT fint F77_FUNC(ipopenoutputfile, IPOPENOUTPUTFILE)(
   fptr* FProblem,
   char* FILENAME,
   fint* PRINTLEVEL,
   int   flen
)
{
   char* filename;
   FUserData* fuser_data = (FUserData*) *FProblem;
   Int printlevel = *PRINTLEVEL;
   fint retval;

   filename = f2cstr(FILENAME, flen);

   retval = OpenIpoptOutputFile(fuser_data->Problem, filename, printlevel);

   free(filename);

   return retval ? OKRetVal : NotOKRetVal;
}

IPOPTLIB_EXPORT void F77_FUNC(ipsetcallback, IPSETCALLBACK)(
   fptr*         FProblem,
   FIntermediate_CB inter_cb
)
{
   FUserData* fuser_data = (FUserData*) *FProblem;
   fuser_data->INTERMEDIATE_CB = inter_cb;
   SetIntermediateCallback(fuser_data->Problem, intermediate_cb);
}

IPOPTLIB_EXPORT void F77_FUNC(ipunsetcallback, IPUNSETCALLBACK)(
   fptr* FProblem
)
{
   FUserData* fuser_data = (FUserData*) *FProblem;
   fuser_data->INTERMEDIATE_CB = NULL;
   SetIntermediateCallback(fuser_data->Problem, NULL);
}
