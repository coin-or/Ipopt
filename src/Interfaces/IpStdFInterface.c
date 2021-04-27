/* Copyright (C) 2004, 2010 International Business Machines and others.
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * Authors:  Carl Laird, Andreas Waechter     IBM    2004-09-03
 */

#include "IpStdCInterface.h"
#include "IpoptConfig.h"
#include "IpTypes.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* if we build without Fortran compiler around, then make up a Fortran name mangling scheme that often works */
#ifndef F77_FUNC
#define F77_FUNC(name,NAME) name ## _
#endif

/* in configure, we checked whether an int* is 32 or 64bit long to decide how much space
 * is needed to store a pointer
 * thus, we can use a void* here to represent a pointer for Fortran
 */
typedef void* fptr;

/** Return value for indicating that evaluation could be done without problem. */
static const ipfint OKRetVal = 0;
/** Return value for indicating that evaluation could not be done without problem. */
static const ipfint NotOKRetVal = 1;

/* Function pointer types for the Fortran callback functions */
typedef void (*FEval_F_CB)(
   ipfint*    N,
   ipnumber*  X,
   ipfint*    NEW_X,
   ipnumber*  OBJVAL,
   ipfint*    IDAT,
   ipnumber*  DDAT,
   ipfint*    IERR
);

typedef void (*FEval_G_CB)(
   ipfint*    N,
   ipnumber*  X,
   ipfint*    NEW_X,
   ipfint*    M,
   ipnumber*  G,
   ipfint*    IDAT,
   ipnumber*  DDAT,
   ipfint*    IERR
);

typedef void (*FEval_Grad_F_CB)(
   ipfint*    N,
   ipnumber*  X,
   ipfint*    NEW_X,
   ipnumber*  GRAD,
   ipfint*    IDAT,
   ipnumber*  DDAT,
   ipfint*    IERR
);

typedef void (*FEval_Jac_G_CB)(
   ipfint*    TASK,
   ipfint*    N,
   ipnumber*  X,
   ipfint*    NEW_X,
   ipfint*    M,
   ipfint*    NNZJAC,
   ipfint*    IROW,
   ipfint*    JCOL,
   ipnumber*  VALUES,
   ipfint*    IDAT,
   ipnumber*  DDAT,
   ipfint*    IERR
);

typedef void (*FEval_Hess_CB)(
   ipfint*    TASK,
   ipfint*    N,
   ipnumber*  X,
   ipfint*    NEW_X,
   ipnumber*  OBJFACT,
   ipfint*    M,
   ipnumber*  LAMBDA,
   ipfint*    NEW_LAM,
   ipfint*    NNZHESS,
   ipfint*    IROW,
   ipfint*    JCOL,
   ipnumber*  VALUES,
   ipfint*    IDAT,
   ipnumber*  DDAT,
   ipfint*    IERR
);

typedef void (*FIntermediate_CB)(
   ipfint*    ALG_MODE,
   ipfint*    ITER_COUNT,
   ipnumber*  OBJVAL,
   ipnumber*  INF_PR,
   ipnumber*  INF_DU,
   ipnumber*  MU,
   ipnumber*  DNORM,
   ipnumber*  REGU_SIZE,
   ipnumber*  ALPHA_DU,
   ipnumber*  ALPHA_PR,
   ipfint*    LS_TRIAL,
   ipfint*    IDAT,
   ipnumber*  DDAT,
   ipfint*    ISTOP
);

struct _FUserData
{
   ipfint*          IDAT;
   ipnumber*        DDAT;
   FEval_F_CB       EVAL_F;
   FEval_G_CB       EVAL_G;
   FEval_Grad_F_CB  EVAL_GRAD_F;
   FEval_Jac_G_CB   EVAL_JAC_G;
   FEval_Hess_CB    EVAL_HESS;
   FIntermediate_CB INTERMEDIATE_CB;
   IpoptProblem     Problem;
};
typedef struct _FUserData FUserData;

static Bool eval_f(
   Index       n,
   Number*     x,
   Bool        new_x,
   Number*     obj_value,
   UserDataPtr user_data
)
{
   ipfint N = n;
   ipfint NEW_X = new_x;
   FUserData* fuser_data = (FUserData*) user_data;
   ipfint* IDAT = fuser_data->IDAT;
   ipnumber*  DDAT = fuser_data->DDAT;
   ipfint IERR = 0;

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
   ipfint N = n;
   ipfint NEW_X = new_x;
   FUserData* fuser_data = (FUserData*) user_data;
   ipfint* IDAT = fuser_data->IDAT;
   ipnumber*  DDAT = fuser_data->DDAT;
   ipfint IERR = 0;

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
   ipfint N = n;
   ipfint NEW_X = new_x;
   ipfint M = m;
   FUserData* fuser_data = (FUserData*) user_data;
   ipfint* IDAT = fuser_data->IDAT;
   ipnumber* DDAT = fuser_data->DDAT;
   ipfint IERR = 0;

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
   ipfint N = n;
   ipfint NEW_X = new_x;
   ipfint M = m;
   ipfint NNZJAC = nele_jac;
   ipfint TASK;
   FUserData* fuser_data = (FUserData*) user_data;
   ipfint* IDAT = fuser_data->IDAT;
   ipnumber*  DDAT = fuser_data->DDAT;
   ipfint IERR = 0;

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
   ipfint N = n;
   ipfint NEW_X = new_x;
   ipfint M = m;
   ipfint NEW_LAM = new_lambda;
   ipfint NNZHESS = nele_hess;
   ipfint TASK;
   FUserData* fuser_data = (FUserData*) user_data;
   ipfint* IDAT = fuser_data->IDAT;
   ipnumber* DDAT = fuser_data->DDAT;
   ipfint IERR = 0;

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
   ipfint ALG_MODE = alg_mod;
   ipfint ITER_COUNT = iter_count;
   ipnumber OBJVAL = obj_value;
   ipnumber INF_PR = inf_pr;
   ipnumber INF_DU = inf_du;
   ipnumber MU = mu;
   ipnumber DNORM = d_norm;
   ipnumber REGU_SIZE = regularization_size;
   ipnumber ALPHA_DU = alpha_du;
   ipnumber ALPHA_PR = alpha_pr;
   ipfint LS_TRIAL = ls_trials;
   ipfint* IDAT = fuser_data->IDAT;
   ipnumber* DDAT = fuser_data->DDAT;
   ipfint ISTOP = 0;

   if( !fuser_data->INTERMEDIATE_CB )
   {
      return (Bool) TRUE;
   }

   fuser_data->INTERMEDIATE_CB(&ALG_MODE, &ITER_COUNT, &OBJVAL, &INF_PR, &INF_DU, &MU, &DNORM, &REGU_SIZE, &ALPHA_DU,
                               &ALPHA_PR, &LS_TRIAL, IDAT, DDAT, &ISTOP);

   return (Bool) (ISTOP == OKRetVal);
}

IPOPTLIB_EXPORT fptr F77_FUNC(ipcreate, IPCREATE)(
   ipfint*         N,
   ipnumber*       X_L,
   ipnumber*       X_U,
   ipfint*         M,
   ipnumber*       G_L,
   ipnumber*       G_U,
   ipfint*         NELE_JAC,
   ipfint*         NELE_HESS,
   ipfint*         IDX_STY,
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

IPOPTLIB_EXPORT ipfint F77_FUNC(ipsolve, IPSOLVE)(
   fptr*      FProblem,
   ipnumber*  X,
   ipnumber*  G,
   ipnumber*  OBJ_VAL,
   ipnumber*  MULT_G,
   ipnumber*  MULT_X_L,
   ipnumber*  MULT_X_U,
   ipfint*    IDAT,
   ipnumber*  DDAT
)
{
   FUserData* fuser_data = (FUserData*) *FProblem;
   UserDataPtr user_data;

   fuser_data->IDAT = IDAT;
   fuser_data->DDAT = DDAT;
   user_data = (UserDataPtr) fuser_data;

   return (ipfint)IpoptSolve(fuser_data->Problem, X, G, OBJ_VAL, MULT_G, MULT_X_L, MULT_X_U, user_data);
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
IPOPTLIB_EXPORT ipfint F77_FUNC(ipaddstroption, IPADDSTROPTION)(
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
   ipfint retval;

   keyword = f2cstr(KEYWORD, klen);
   val = f2cstr(VALUE, vlen);

   retval = AddIpoptStrOption(fuser_data->Problem, keyword, val);

   free(val);
   free(keyword);

   return retval ? OKRetVal : NotOKRetVal;
}

IPOPTLIB_EXPORT ipfint F77_FUNC(ipaddnumoption, IPADDNUMOPTION)(
   fptr*     FProblem,
   char*     KEYWORD,
   ipnumber* VALUE,
   int       klen
)
{
   char* keyword;
   FUserData* fuser_data = (FUserData*) *FProblem;
   ipfint retval;

   keyword = f2cstr(KEYWORD, klen);

   retval = AddIpoptNumOption(fuser_data->Problem, keyword, *VALUE);

   free(keyword);

   return retval ? OKRetVal : NotOKRetVal;
}

IPOPTLIB_EXPORT ipfint F77_FUNC(ipaddintoption, IPADDINTOPTION)(
   fptr*   FProblem,
   char*   KEYWORD,
   ipfint* VALUE,
   int     klen
)
{
   char* keyword;
   FUserData* fuser_data = (FUserData*) *FProblem;
   Int value = *VALUE;
   ipfint retval;

   keyword = f2cstr(KEYWORD, klen);

   retval = AddIpoptIntOption(fuser_data->Problem, keyword, value);

   free(keyword);

   return retval ? OKRetVal : NotOKRetVal;
}

IPOPTLIB_EXPORT ipfint F77_FUNC(ipopenoutputfile, IPOPENOUTPUTFILE)(
   fptr*   FProblem,
   char*   FILENAME,
   ipfint* PRINTLEVEL,
   int     flen
)
{
   char* filename;
   FUserData* fuser_data = (FUserData*) *FProblem;
   Int printlevel = *PRINTLEVEL;
   ipfint retval;

   filename = f2cstr(FILENAME, flen);

   retval = OpenIpoptOutputFile(fuser_data->Problem, filename, printlevel);

   free(filename);

   return retval ? OKRetVal : NotOKRetVal;
}

IPOPTLIB_EXPORT void F77_FUNC(ipsetcallback, IPSETCALLBACK)(
   fptr*            FProblem,
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

IPOPTLIB_EXPORT ipfint F77_FUNC(ipgetcurriterate, IPGETCURRITERATE)(
   fptr*      FProblem,
   ipfint*    scaled,
   ipfint*    get_X,
   ipfint*    get_Z,
   ipfint*    get_G,
   ipfint*    get_LAMBDA,
   ipfint*    n,
   ipnumber*  X,
   ipnumber*  Z_L,
   ipnumber*  Z_U,
   ipfint*    m,
   ipnumber*  G,
   ipnumber*  LAMBDA
)
{
   FUserData* fuser_data = (FUserData*) *FProblem;
   return GetIpoptCurrentIterate(fuser_data->Problem, *scaled != 0,
      *n, *get_X ? X : NULL, *get_Z ? Z_L : NULL, *get_Z ? Z_U : NULL,
      *m, *get_G ? G : NULL, *get_LAMBDA ? LAMBDA : NULL) ? OKRetVal : NotOKRetVal;
}

IPOPTLIB_EXPORT ipfint F77_FUNC(ipgetcurrviolations, IPGETCURRVIOLATIONS)(
   fptr*      FProblem,
   ipfint*    scaled,
   ipfint*    get_bound_violation,
   ipfint*    get_compl,
   ipfint*    get_grad_lag_x,
   ipfint*    get_nlp_constraint_violation,
   ipfint*    n,
   ipnumber*  X_L_VIOLATION,
   ipnumber*  X_U_VIOLATION,
   ipnumber*  COMPL_X_L,
   ipnumber*  COMPL_X_U,
   ipnumber*  GRAD_LAG_X,
   ipfint*    m,
   ipnumber*  NLP_CONSTRAINT_VIOLATION,
   ipnumber*  COMPL_G
)
{
   FUserData* fuser_data = (FUserData*) *FProblem;
   return GetIpoptCurrentViolations(fuser_data->Problem, *scaled != 0,
      *n,
      *get_bound_violation ? X_L_VIOLATION : NULL, *get_bound_violation ? X_U_VIOLATION : NULL,
      *get_compl ? COMPL_X_L : NULL, *get_compl ? COMPL_X_U : NULL,
      *get_grad_lag_x ? GRAD_LAG_X : NULL,
      *m,
      *get_nlp_constraint_violation ? NLP_CONSTRAINT_VIOLATION : NULL,
      *get_compl ? COMPL_G : NULL)
      ? OKRetVal : NotOKRetVal;
}
