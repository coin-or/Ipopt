/********************************************************************
   Copyright (C) 2004, International Business Machines and others.
   All Rights Reserved.
   This code is published under the Common Public License.

   $Id$

   Authors:  Carl Laird, Andreas Waechter     IBM    2004-09-03
 ********************************************************************/

#include "IpStdCInterface.h"
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* ToDo: The following needs to be adapted based on configuration */
typedef int fint;
typedef double fdouble;
typedef long fptr;

/* Function pointer types for the Fortran callback functions */
typedef void (*FEval_F_CB)(fint* N, fdouble* X, fint* NEW_X,
			   fdouble* OBJVAL, fint* IDAT, fdouble* DDAT,
			   fint* IERR);
typedef void (*FEval_G_CB)(fint* N, fdouble* X, fint* NEW_X,
			   fint* M, fdouble* G, fint* IDAT, fdouble* DDAT,
			   fint* IERR);
typedef void (*FEval_Grad_F_CB)(fint *N, fdouble* X, fint* NEW_X,
				fdouble* GRAD, fint* IDAT, fdouble* DDAT,
				fint* IERR);
typedef void (*FEval_Jac_G_CB)(fint* TASK, fint* N, fdouble* X, fint* NEW_X,
			       fint* M, fint* NNZJAC, fint* IROW, fint* JCOL,
			       fdouble* VALUES, fint* IDAT, fdouble* DDAT,
			       fint* IERR);
typedef void (*FEval_Hess_CB)(fint* TASK, fint* N, fdouble* X, fint* NEW_X,
			      fdouble *OBJFACT, fint* M, fdouble* LAMBDA,
			      fint* NEW_LAM, fint* NNZHESS, fint* IROW,
			      fint* JCOL, fdouble* VALUES, fint* IDAT,
			      fdouble* DDAT, fint* IERR);

static const fint OKRetVal = 0;

struct _FUserData {
  fint* IDAT;
  fdouble* DDAT;
  FEval_F_CB EVAL_F;
  FEval_G_CB EVAL_G;
  FEval_Grad_F_CB EVAL_GRAD_F;
  FEval_Jac_G_CB EVAL_JAC_G;
  FEval_Hess_CB EVAL_HESS;
};

typedef struct _FUserData FUserData;

static Bool eval_f(Index n, Number* x, Bool new_x,
		   Number* obj_value, UserDataPtr user_data)
{
  fint N = n;
  fint NEW_X = new_x;
  FUserData* fuser_data = (FUserData*)user_data;
  fint* IDAT = fuser_data->IDAT;
  fdouble* DDAT = fuser_data->DDAT;
  fint IERR=0;

  fuser_data->EVAL_F(&N, x, &NEW_X, obj_value, IDAT, DDAT, &IERR);

  return (Bool) (IERR==OKRetVal);
}

static Bool eval_grad_f(Index n, Number* x, Bool new_x,
			Number* grad_f, UserDataPtr user_data)
{
  fint N = n;
  fint NEW_X = new_x;
  FUserData* fuser_data = (FUserData*)user_data;
  fint* IDAT = fuser_data->IDAT;
  fdouble* DDAT = fuser_data->DDAT;
  fint IERR=0;

  fuser_data->EVAL_GRAD_F(&N, x, &NEW_X, grad_f, IDAT, DDAT, &IERR);

  return (Bool) (IERR==OKRetVal);
}

static Bool eval_g(Index n, Number* x, Bool new_x,
		   Index m, Number* g, UserDataPtr user_data)
{
  fint N = n;
  fint NEW_X = new_x;
  fint M = m;
  FUserData* fuser_data = (FUserData*)user_data;
  fint* IDAT = fuser_data->IDAT;
  fdouble* DDAT = fuser_data->DDAT;
  fint IERR=0;

  fuser_data->EVAL_G(&N, x, &NEW_X, &M, g, IDAT, DDAT, &IERR);

  return (Bool) (IERR==OKRetVal);
}

static Bool eval_jac_g(Index n, Number *x, Bool new_x, Index m, Index nele_jac,
		       Index *iRow, Index *jCol, Number *values,
		       UserDataPtr user_data)
{
  fint N = n;
  fint NEW_X = new_x;
  fint M = m;
  fint NNZJAC = nele_jac;
  fint TASK;
  FUserData* fuser_data = (FUserData*)user_data;
  fint* IDAT = fuser_data->IDAT;
  fdouble* DDAT = fuser_data->DDAT;
  fint IERR=0;

  if(iRow && jCol && !values) {
    /* Only request the structure */
    TASK = 0;
  }
  else if (!iRow && !jCol && values) {
    /* Only request the values */
    TASK = 1;
  } else {
    printf("Error in IpStdFInterface eval_jac_g!\n");
    return (Bool) 0;
  }

  fuser_data->EVAL_JAC_G(&TASK, &N, x, &NEW_X, &M, &NNZJAC, iRow, jCol,
			 values, IDAT, DDAT, &IERR);

  return (Bool) (IERR==OKRetVal);
}

static Bool eval_h(Index n, Number *x, Bool new_x, Number obj_factor,
		   Index m, Number *lambda, Bool new_lambda, 
		   Index nele_hess, Index *iRow, Index *jCol,
		   Number *values, UserDataPtr user_data)
{
  fint N = n;
  fint NEW_X = new_x;
  fint M = m;
  fint NEW_LAM = new_lambda;
  fint NNZHESS = nele_hess;
  fint TASK;
  FUserData* fuser_data = (FUserData*)user_data;
  fint* IDAT = fuser_data->IDAT;
  fdouble* DDAT = fuser_data->DDAT;
  fint IERR=0;

  if(iRow && jCol && !values) {
    /* Only request the structure */
    TASK = 0;
  }
  else if (!iRow && !jCol && values) {
    /* Only request the values */
    TASK = 1;
  } else {
    printf("Error in IpStdFInterface eval_hess!\n");
    return (Bool) 0;
  }

  fuser_data->EVAL_HESS(&TASK, &N, x, &NEW_X, &obj_factor,
			&M, lambda, &NEW_LAM, &NNZHESS, iRow, jCol,
			values, IDAT, DDAT, &IERR);

  return (Bool) (IERR==OKRetVal);
}

void F77_FUNC(ipopt,IPOPT)(fint* N,
			   fdouble* X,
			   fdouble* X_L,
			   fdouble* X_U,
			   fint* M,
			   fdouble* G,
			   fdouble* G_L,
			   fdouble* G_U,
			   fint* NELE_JAC,
			   fint* NELE_HESS,
			   fdouble* OBJ_VAL,
			   fdouble* MULT_G,
			   fdouble* MULT_X_L,
			   fdouble* MULT_X_U,
			   FEval_F_CB EVAL_F,
			   FEval_G_CB EVAL_G,
			   FEval_Grad_F_CB EVAL_GRAD_F,
			   FEval_Jac_G_CB EVAL_JAC_G,
			   FEval_Hess_CB EVAL_HESS,
			   fptr* OPTIONS,
			   fint* IDAT,
			   fdouble* DDAT,
			   fint* IERR)
{
  Index n = *N;
  Index m = *M;
  Index nele_jac = *NELE_JAC;
  Index nele_hess = *NELE_HESS;
  OptionsPtr options = (OptionsPtr)*OPTIONS;

  FUserData* fuser_data;
  UserDataPtr user_data;

  fuser_data = (FUserData*) malloc(sizeof(FUserData));
  fuser_data->IDAT = IDAT;
  fuser_data->DDAT = DDAT;
  fuser_data->EVAL_F = EVAL_F;
  fuser_data->EVAL_G = EVAL_G;
  fuser_data->EVAL_GRAD_F = EVAL_GRAD_F;
  fuser_data->EVAL_JAC_G = EVAL_JAC_G;
  fuser_data->EVAL_HESS = EVAL_HESS;
  user_data = (UserDataPtr) fuser_data;

  *IERR = (fint)IpoptSolve(n, X, X_L, X_U, m, G, G_L, G_U, nele_jac,
			   nele_hess, OBJ_VAL, MULT_G, MULT_X_L, MULT_X_U,
			   eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h,
			   options, user_data);

  free(fuser_data);

  return;
}

fptr F77_FUNC(ipnewopts,IPNEWOPTS)()
{
  return (fptr) Ipopt_NewOptions();
}

void F77_FUNC(ipdelopts,IPDELOPTS)(fptr* OPTIONS)
{
  if ((OptionsPtr)*OPTIONS!=NULL) {
    Ipopt_DeleteOptions((OptionsPtr)*OPTIONS);
  }
}

static char* f2cstr(char* FSTR, int slen) {
  int len;
  char* cstr;
  for (len=slen;len>0;len--) {
    if (FSTR[len-1]!=' ') {
      break;
    }
  }
  cstr = (char*)malloc(sizeof(char)*(len+1));
  strncpy(cstr, FSTR, len);
  cstr[len+1]='\0';

  return cstr;
}

void F77_FUNC(ipaddopt,IPADDOPT)(fptr* OPTIONS, char* VALUE,
				 char* KEYWORD, int vlen, int klen)
{
  char* keyword;
  char* val;

  keyword = f2cstr(KEYWORD, klen);
  val = f2cstr(VALUE, vlen);

  Ipopt_AddOption((OptionsPtr)*OPTIONS, keyword, val);

  free(val);
  free(keyword);
}

void F77_FUNC(ipaddnumopt,IPADDNUMOPT)(fptr* OPTIONS, fdouble* VAL,
				       char* KEYWORD, int klen)
{
  char* keyword;

  keyword = f2cstr(KEYWORD, klen);

  Ipopt_AddNumOption((OptionsPtr)*OPTIONS, keyword, *VAL);

  free(keyword);
}

void F77_FUNC(ipaddintopt,IPADDINTOPT)(fptr* OPTIONS, fint* VAL,
				       char* KEYWORD, int klen)
{
  char* keyword;

  keyword = f2cstr(KEYWORD, klen);

  Ipopt_AddIntOption((OptionsPtr)*OPTIONS, keyword, *VAL);

  free(keyword);
}
