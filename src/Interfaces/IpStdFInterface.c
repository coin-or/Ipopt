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
static const ipindex OKRetVal = 0;
/** Return value for indicating that evaluation could not be done without problem. */
static const ipindex NotOKRetVal = 1;

/* Function pointer types for the Fortran callback functions */
typedef void (*FEval_F_CB)(
   ipindex*    N,
   ipnumber*   X,
   ipindex*    NEW_X,
   ipnumber*   OBJVAL,
   ipindex*    IDAT,
   ipnumber*   DDAT,
   ipindex*    IERR
);

typedef void (*FEval_G_CB)(
   ipindex*    N,
   ipnumber*   X,
   ipindex*    NEW_X,
   ipindex*    M,
   ipnumber*   G,
   ipindex*    IDAT,
   ipnumber*   DDAT,
   ipindex*    IERR
);

typedef void (*FEval_Grad_F_CB)(
   ipindex*    N,
   ipnumber*   X,
   ipindex*    NEW_X,
   ipnumber*   GRAD,
   ipindex*    IDAT,
   ipnumber*   DDAT,
   ipindex*    IERR
);

typedef void (*FEval_Jac_G_CB)(
   ipindex*    TASK,
   ipindex*    N,
   ipnumber*   X,
   ipindex*    NEW_X,
   ipindex*    M,
   ipindex*    NNZJAC,
   ipindex*    IROW,
   ipindex*    JCOL,
   ipnumber*   VALUES,
   ipindex*    IDAT,
   ipnumber*   DDAT,
   ipindex*    IERR
);

typedef void (*FEval_Hess_CB)(
   ipindex*    TASK,
   ipindex*    N,
   ipnumber*   X,
   ipindex*    NEW_X,
   ipnumber*   OBJFACT,
   ipindex*    M,
   ipnumber*   LAMBDA,
   ipindex*    NEW_LAM,
   ipindex*    NNZHESS,
   ipindex*    IROW,
   ipindex*    JCOL,
   ipnumber*   VALUES,
   ipindex*    IDAT,
   ipnumber*   DDAT,
   ipindex*    IERR
);

typedef void (*FIntermediate_CB)(
   ipindex*    ALG_MODE,
   ipindex*    ITER_COUNT,
   ipnumber*   OBJVAL,
   ipnumber*   INF_PR,
   ipnumber*   INF_DU,
   ipnumber*   MU,
   ipnumber*   DNORM,
   ipnumber*   REGU_SIZE,
   ipnumber*   ALPHA_DU,
   ipnumber*   ALPHA_PR,
   ipindex*    LS_TRIAL,
   ipindex*    IDAT,
   ipnumber*   DDAT,
   ipindex*    ISTOP
);

struct _FUserData
{
   ipindex*         IDAT;
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

static bool eval_f(
   ipindex     n,
   ipnumber*   x,
   bool        new_x,
   ipnumber*   obj_value,
   UserDataPtr user_data
)
{
   ipindex N = n;
   ipindex NEW_X = new_x;
   FUserData* fuser_data = (FUserData*) user_data;
   ipindex* IDAT = fuser_data->IDAT;
   ipnumber*  DDAT = fuser_data->DDAT;
   ipindex IERR = 0;

   fuser_data->EVAL_F(&N, x, &NEW_X, obj_value, IDAT, DDAT, &IERR);

   return IERR == OKRetVal;
}

static bool eval_grad_f(
   ipindex     n,
   ipnumber*   x,
   bool        new_x,
   ipnumber*   grad_f,
   UserDataPtr user_data
)
{
   ipindex N = n;
   ipindex NEW_X = new_x;
   FUserData* fuser_data = (FUserData*) user_data;
   ipindex* IDAT = fuser_data->IDAT;
   ipnumber*  DDAT = fuser_data->DDAT;
   ipindex IERR = 0;

   fuser_data->EVAL_GRAD_F(&N, x, &NEW_X, grad_f, IDAT, DDAT, &IERR);

   return IERR == OKRetVal;
}

static bool eval_g(
   ipindex     n,
   ipnumber*   x,
   bool        new_x,
   ipindex     m,
   ipnumber*   g,
   UserDataPtr user_data
)
{
   ipindex N = n;
   ipindex NEW_X = new_x;
   ipindex M = m;
   FUserData* fuser_data = (FUserData*) user_data;
   ipindex* IDAT = fuser_data->IDAT;
   ipnumber* DDAT = fuser_data->DDAT;
   ipindex IERR = 0;

   fuser_data->EVAL_G(&N, x, &NEW_X, &M, g, IDAT, DDAT, &IERR);

   return IERR == OKRetVal;
}

static bool eval_jac_g(
   ipindex     n,
   ipnumber*   x,
   bool        new_x,
   ipindex     m,
   ipindex     nele_jac,
   ipindex*    iRow,
   ipindex*    jCol,
   ipnumber*   values,
   UserDataPtr user_data
)
{
   ipindex N = n;
   ipindex NEW_X = new_x;
   ipindex M = m;
   ipindex NNZJAC = nele_jac;
   ipindex TASK;
   FUserData* fuser_data = (FUserData*) user_data;
   ipindex* IDAT = fuser_data->IDAT;
   ipnumber*  DDAT = fuser_data->DDAT;
   ipindex IERR = 0;

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
      return false;
   }

   fuser_data->EVAL_JAC_G(&TASK, &N, x, &NEW_X, &M, &NNZJAC, iRow, jCol, values, IDAT, DDAT, &IERR);

   return IERR == OKRetVal;
}

static bool eval_h(
   ipindex     n,
   ipnumber*   x,
   bool        new_x,
   ipnumber    obj_factor,
   ipindex     m,
   ipnumber*   lambda,
   bool        new_lambda,
   ipindex     nele_hess,
   ipindex*    iRow,
   ipindex*    jCol,
   ipnumber*   values,
   UserDataPtr user_data
)
{
   ipindex N = n;
   ipindex NEW_X = new_x;
   ipindex M = m;
   ipindex NEW_LAM = new_lambda;
   ipindex NNZHESS = nele_hess;
   ipindex TASK;
   FUserData* fuser_data = (FUserData*) user_data;
   ipindex* IDAT = fuser_data->IDAT;
   ipnumber* DDAT = fuser_data->DDAT;
   ipindex IERR = 0;

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
      return false;
   }

   fuser_data->EVAL_HESS(&TASK, &N, x, &NEW_X, &obj_factor, &M, lambda, &NEW_LAM, &NNZHESS, iRow, jCol, values, IDAT, DDAT, &IERR);

   return IERR == OKRetVal;
}

static bool intermediate_cb(
   ipindex     alg_mod,
   ipindex     iter_count,
   ipnumber    obj_value,
   ipnumber    inf_pr,
   ipnumber    inf_du,
   ipnumber    mu,
   ipnumber    d_norm,
   ipnumber    regularization_size,
   ipnumber    alpha_du,
   ipnumber    alpha_pr,
   ipindex     ls_trials,
   UserDataPtr user_data
)
{
   FUserData* fuser_data = (FUserData*) user_data;
   ipindex ALG_MODE = alg_mod;
   ipindex ITER_COUNT = iter_count;
   ipnumber OBJVAL = obj_value;
   ipnumber INF_PR = inf_pr;
   ipnumber INF_DU = inf_du;
   ipnumber MU = mu;
   ipnumber DNORM = d_norm;
   ipnumber REGU_SIZE = regularization_size;
   ipnumber ALPHA_DU = alpha_du;
   ipnumber ALPHA_PR = alpha_pr;
   ipindex LS_TRIAL = ls_trials;
   ipindex* IDAT = fuser_data->IDAT;
   ipnumber* DDAT = fuser_data->DDAT;
   ipindex ISTOP = 0;

   if( !fuser_data->INTERMEDIATE_CB )
   {
      return true;
   }

   fuser_data->INTERMEDIATE_CB(&ALG_MODE, &ITER_COUNT, &OBJVAL, &INF_PR, &INF_DU, &MU, &DNORM, &REGU_SIZE, &ALPHA_DU,
                               &ALPHA_PR, &LS_TRIAL, IDAT, DDAT, &ISTOP);

   return ISTOP == OKRetVal;
}

IPOPTLIB_EXPORT fptr F77_FUNC(ipcreate, IPCREATE)(
   ipindex*        N,
   ipnumber*       X_L,
   ipnumber*       X_U,
   ipindex*        M,
   ipnumber*       G_L,
   ipnumber*       G_U,
   ipindex*        NELE_JAC,
   ipindex*        NELE_HESS,
   ipindex*        IDX_STY,
   FEval_F_CB      EVAL_F,
   FEval_G_CB      EVAL_G,
   FEval_Grad_F_CB EVAL_GRAD_F,
   FEval_Jac_G_CB  EVAL_JAC_G,
   FEval_Hess_CB   EVAL_HESS
)
{
   ipindex n = *N;
   ipindex m = *M;
   ipindex nele_jac = *NELE_JAC;
   ipindex nele_hess = *NELE_HESS;
   ipindex index_style = *IDX_STY;

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

IPOPTLIB_EXPORT ipindex F77_FUNC(ipsolve, IPSOLVE)(
   fptr*      FProblem,
   ipnumber*  X,
   ipnumber*  G,
   ipnumber*  OBJ_VAL,
   ipnumber*  MULT_G,
   ipnumber*  MULT_X_L,
   ipnumber*  MULT_X_U,
   ipindex*   IDAT,
   ipnumber*  DDAT
)
{
   FUserData* fuser_data = (FUserData*) *FProblem;
   UserDataPtr user_data;

   fuser_data->IDAT = IDAT;
   fuser_data->DDAT = DDAT;
   user_data = (UserDataPtr) fuser_data;

   return IpoptSolve(fuser_data->Problem, X, G, OBJ_VAL, MULT_G, MULT_X_L, MULT_X_U, user_data);
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
   if( cstr != NULL )
   {
      strncpy(cstr, FSTR, len);
      cstr[len] = '\0';
   }

   return cstr;
}

/* ToDo make sure position of vlen and klen are at the right place */
IPOPTLIB_EXPORT ipindex F77_FUNC(ipaddstroption, IPADDSTROPTION)(
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
   ipindex retval;

   keyword = f2cstr(KEYWORD, klen);
   val = f2cstr(VALUE, vlen);

   retval = AddIpoptStrOption(fuser_data->Problem, keyword, val);

   free(val);
   free(keyword);

   return retval ? OKRetVal : NotOKRetVal;
}

IPOPTLIB_EXPORT ipindex F77_FUNC(ipaddnumoption, IPADDNUMOPTION)(
   fptr*     FProblem,
   char*     KEYWORD,
   ipnumber* VALUE,
   int       klen
)
{
   char* keyword;
   FUserData* fuser_data = (FUserData*) *FProblem;
   ipindex retval;

   keyword = f2cstr(KEYWORD, klen);

   retval = AddIpoptNumOption(fuser_data->Problem, keyword, *VALUE);

   free(keyword);

   return retval ? OKRetVal : NotOKRetVal;
}

IPOPTLIB_EXPORT ipindex F77_FUNC(ipaddintoption, IPADDINTOPTION)(
   fptr*   FProblem,
   char*   KEYWORD,
   ipindex* VALUE,
   int     klen
)
{
   char* keyword;
   FUserData* fuser_data = (FUserData*) *FProblem;
   ipindex value = *VALUE;
   ipindex retval;

   keyword = f2cstr(KEYWORD, klen);

   retval = AddIpoptIntOption(fuser_data->Problem, keyword, value);

   free(keyword);

   return retval ? OKRetVal : NotOKRetVal;
}

IPOPTLIB_EXPORT ipindex F77_FUNC(ipopenoutputfile, IPOPENOUTPUTFILE)(
   fptr*    FProblem,
   char*    FILENAME,
   ipindex* PRINTLEVEL,
   int      flen
)
{
   char* filename;
   FUserData* fuser_data = (FUserData*) *FProblem;
   int printlevel = *PRINTLEVEL;
   ipindex retval;

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

/// @since 3.14.0
IPOPTLIB_EXPORT ipindex F77_FUNC(ipgetcurriterate, IPGETCURRITERATE)(
   fptr*      FProblem,
   ipindex*   scaled,
   ipindex*   get_X,
   ipindex*   get_Z,
   ipindex*   get_G,
   ipindex*   get_LAMBDA,
   ipindex*   n,
   ipnumber*  X,
   ipnumber*  Z_L,
   ipnumber*  Z_U,
   ipindex*   m,
   ipnumber*  G,
   ipnumber*  LAMBDA
)
{
   FUserData* fuser_data = (FUserData*) *FProblem;
   return GetIpoptCurrentIterate(fuser_data->Problem, *scaled != 0,
                                 *n, *get_X ? X : NULL, *get_Z ? Z_L : NULL, *get_Z ? Z_U : NULL,
                                 *m, *get_G ? G : NULL, *get_LAMBDA ? LAMBDA : NULL) ? OKRetVal : NotOKRetVal;
}

/// @since 3.14.0
IPOPTLIB_EXPORT ipindex F77_FUNC(ipgetcurrviolations, IPGETCURRVIOLATIONS)(
   fptr*       FProblem,
   ipindex*    scaled,
   ipindex*    get_bound_violation,
   ipindex*    get_compl,
   ipindex*    get_grad_lag_x,
   ipindex*    get_nlp_constraint_violation,
   ipindex*    n,
   ipnumber*   X_L_VIOLATION,
   ipnumber*   X_U_VIOLATION,
   ipnumber*   COMPL_X_L,
   ipnumber*   COMPL_X_U,
   ipnumber*   GRAD_LAG_X,
   ipindex*    m,
   ipnumber*   NLP_CONSTRAINT_VIOLATION,
   ipnumber*   COMPL_G
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

/// @since 3.14.7
IPOPTLIB_EXPORT ipindex F77_FUNC(ipsetproblemscaling, IPSETPROBLEMSCALING)(
   fptr*      FProblem,
   ipnumber*  obj_scaling,
   ipnumber*  X_SCALING,
   ipnumber*  G_SCALING
)
{
   FUserData* fuser_data = (FUserData*) *FProblem;
   return SetIpoptProblemScaling(fuser_data->Problem, *obj_scaling, X_SCALING, G_SCALING) ? OKRetVal : NotOKRetVal;
}
