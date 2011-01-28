/* Copyright (C) 2009 International Business Machines.
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * $Id$
 *
 * Author:  Andreas Waechter               IBM    2009-04-02
 */

/*
// This file is part of the Ipopt tutorial.  It is a version with
// mistakes for the C implemention of the coding exercise problem (in
// AMPL formulation):
//
// param n := 4;
//
// var x {1..n} <= 0, >= -1.5, := -0.5;
//
// minimize obj:
//   sum{i in 1..n} (x[i]-1)^2;
//
// subject to constr {i in 2..n-1}:
//   (x[i]^2+1.5*x[i]-i/n)*cos(x[i+1]) - x[i-1] = 0;
//
// The constant term "i/n" in the constraint is supposed to be input data
*/

#include "IpStdCInterface.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>

/* Function Declarations */
Bool eval_f(Index n, Number* x, Bool new_x,
            Number* obj_value, UserDataPtr user_data);

Bool eval_grad_f(Index n, Number* x, Bool new_x,
                 Number* grad_f, UserDataPtr user_data);

Bool eval_g(Index n, Number* x, Bool new_x,
            Index m, Number* g, UserDataPtr user_data);

Bool eval_jac_g(Index n, Number *x, Bool new_x,
                Index m, Index nele_jac,
                Index *iRow, Index *jCol, Number *values,
                UserDataPtr user_data);

Bool eval_h(Index n, Number *x, Bool new_x, Number obj_factor,
            Index m, Number *lambda, Bool new_lambda,
            Index nele_hess, Index *iRow, Index *jCol,
            Number *values, UserDataPtr user_data);

/* Structure to communicate problem data */
struct _ProblemData {
  int N;
  double* a;
};
typedef struct _ProblemData* ProblemData;

/* Main Program */
int main()
{
  Index n=-1;                          /* number of variables */
  Index m=-1;                          /* number of constraints */
  Index nele_jac;                      /* number of nonzeros in Jacobian */
  Index nele_hess;                     /* number of nonzeros in Hessian */
  Index index_style;                   /* indexing style for matrices */
  Number* x_L = NULL;                  /* lower bounds on x */
  Number* x_U = NULL;                  /* upper bounds on x */
  Number* g_L = NULL;                  /* lower bounds on g */
  Number* g_U = NULL;                  /* upper bounds on g */
  IpoptProblem nlp = NULL;             /* IpoptProblem */
  enum ApplicationReturnStatus status; /* Solve return code */
  Number* x = NULL;                    /* starting point and solution vector */
  Number* mult_x_L = NULL;             /* lower bound multipliers
  					  at the solution */
  Number* mult_x_U = NULL;             /* upper bound multipliers
  					  at the solution */
  Number obj;                          /* objective value */
  Index i;                             /* generic counter */

  int size;                            /* Size of the problem */
  ProblemData PD;                      /* Pointer to structure with problem data */

  /* Specify size of the problem */
  size = 5; /* 100; */

  /* Set the problem data */
  PD = (ProblemData)malloc(sizeof(struct _ProblemData));
  PD->N = size;
  PD->a = malloc(sizeof(double)*(size-2));
  for (i=0; i<size-2; i++) {
    PD->a[i] = ((double)(i+2))/(double)size;
  }

  /* set the number of variables and allocate space for the bounds */
  n=size;
  x_L = (Number*)malloc(sizeof(Number)*n);
  x_U = (Number*)malloc(sizeof(Number)*n);
  /* set the values for the variable bounds */
  for (i=0; i<n; i++) {
    x_L[i] = -1.5;
    x_U[i] = 0.;
  }

  /* set the number of constraints and allocate space for the bounds */
  m=size-2;
  g_L = (Number*)malloc(sizeof(Number)*m);
  g_U = (Number*)malloc(sizeof(Number)*m);
  /* set the values of the constraint bounds */
  for (i=0; i<m; i++) {
    g_L[i] = 0.;
    g_U[i] = 0.;
  }
  
  /* Number of nonzeros in the Jacobian of the constraints
     each constraint has three nonzeros */
  nele_jac = 3*m;

  /* Number of nonzeros in the Hessian of the Lagrangian (lower or
     upper triangual part only)
     We have the full diagonal, and the first off-diagonal except for
     the first and last variable */
  nele_hess = n + (n-2);

  /* indexing style for matrices */
  index_style = 0; /* C-style; start counting of rows and column
  			    indices at 0 */

  /* create the IpoptProblem */
  nlp = CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess,
                           index_style, &eval_f, &eval_g, &eval_grad_f,
                           &eval_jac_g, &eval_h);

  /* We can free the memory now - the values for the bounds have been
     copied internally in CreateIpoptProblem */
  free(x_L);
  free(x_U);
  free(g_L);
  free(g_U);

  /* Set some options.  Note the following ones are only examples,
     they might not be suitable for your problem. */
  AddIpoptNumOption(nlp, "tol", 1e-7);
  AddIpoptStrOption(nlp, "mu_strategy", "adaptive");
  AddIpoptStrOption(nlp, "output_file", "ipopt.out");

  /* allocate space for the initial point and set the values */
  x = (Number*)malloc(sizeof(Number)*n);
  for (i=0; i<n; i++) {
    x[i] = -0.5;
  }

  /*#define skip_me */
#ifdef skip_me
  /* If checking derivatives, if is useful to choose different values */
  for (i=0; i<n; i++) {
    x[i] = -0.5+0.1*i/n;
  }
#endif

  /* allocate space to store the bound multipliers at the solution */
  mult_x_L = (Number*)malloc(sizeof(Number)*n);
  mult_x_U = (Number*)malloc(sizeof(Number)*n);

  /* solve the problem */
  status = IpoptSolve(nlp, x, NULL, &obj, NULL, mult_x_L, mult_x_U, (void*)PD);

  if (status == Solve_Succeeded) {
    printf("\n\nSolution of the primal variables, x\n");
    for (i=0; i<n; i++) {
      printf("x[%d] = %e\n", i, x[i]);
    }

    printf("\n\nSolution of the bound multipliers, z_L and z_U\n");
    for (i=0; i<n; i++) {
      printf("z_L[%d] = %e\n", i, mult_x_L[i]);
    }
    for (i=0; i<n; i++) {
      printf("z_U[%d] = %e\n", i, mult_x_U[i]);
    }

    printf("\n\nObjective value\n");
    printf("f(x*) = %e\n", obj);
  }

  /* free allocated memory */
  FreeIpoptProblem(nlp);
  free(x);
  free(mult_x_L);
  free(mult_x_U);

  /* also for our user data */
  free(PD->a);
  free(PD);

  return 0;
}


/* Function Implementations */
Bool eval_f(Index n, Number* x, Bool new_x,
            Number* obj_value, UserDataPtr user_data)
{
  int i;
  ProblemData PD = (ProblemData)user_data;
  assert(n == PD->N);

  *obj_value = 0.;
  for (i=0; i<n; i++) {
    *obj_value += (x[i]-1.)*(x[i]-1.);
  }

  return TRUE;
}

Bool eval_grad_f(Index n, Number* x, Bool new_x,
                 Number* grad_f, UserDataPtr user_data)
{
  int i;
  ProblemData PD = (ProblemData)user_data;
  assert(n == PD->N);

  for (i=1; i<n; i++) {
    grad_f[i] = 2.*(x[i]-1.);
  }

  return TRUE;
}

Bool eval_g(Index n, Number* x, Bool new_x,
            Index m, Number* g, UserDataPtr user_data)
{
  int j;
  ProblemData PD = (ProblemData)user_data;
  double* a = PD->a;
  assert(n == PD->N);
  assert(m == PD->N-2);

  for (j=0; j<m; j++) {
    g[j] = (x[j+1]*x[j+1] + 1.5*x[j+1] -a[j])*cos(x[j+2]) - x[j];
  }

  return TRUE;
}

Bool eval_jac_g(Index n, Number *x, Bool new_x,
                Index m, Index nele_jac,
                Index *iRow, Index *jCol, Number *values,
                UserDataPtr user_data)
{
  int j, inz;
  ProblemData PD = (ProblemData)user_data;
  double* a = PD->a;

  if (values == NULL) {
    /* return the structure of the jacobian */

    inz = 0;
    for (j=0; j<m; j++) {
      iRow[inz] = j;
      jCol[inz] = j;
      inz++;
      iRow[inz] = j;
      jCol[inz] = j+1;
      inz++;
      iRow[inz] = j;
      jCol[inz] = j+1;
      inz++;
    }
    /* sanity check */
    assert(inz==nele_jac);
  }
  else {
    /* return the values of the jacobian of the constraints */

    inz = 0;
    for (j=0; j<m; j++) {
      values[inz] = 1.;
      inz++;
      values[inz] = (2.*x[j+1]+1.5)*cos(x[j+2]);
      inz++;
      values[inz] = -(x[j+1]*x[j+1] + 1.5*x[j+1] -a[j])*sin(x[j+2]);
      inz++;
    }
    /* sanity check */
    assert(inz==nele_jac);
  }

  return TRUE;
}

Bool eval_h(Index n, Number *x, Bool new_x, Number obj_factor,
            Index m, Number *lambda, Bool new_lambda,
            Index nele_hess, Index *iRow, Index *jCol,
            Number *values, UserDataPtr user_data)
{
  int i, inz;
  ProblemData PD = (ProblemData)user_data;
  double* a = PD->a;

  if (values == NULL) {
    /* return the structure. This is a symmetric matrix, fill the
     * upper right triangle only. */

    inz = 0;

    /* First variable has only a diagonal entry */
    iRow[inz] = 0;
    jCol[inz] = 0;
    inz++;

    /* Next ones have first off-diagonal and diagonal */
    for (i=1; i<n-1; i++) {
      iRow[inz] = i;
      jCol[inz] = i;
      inz++;
      iRow[inz] = i;
      jCol[inz] = i+1;
      inz++;
    }

    /* Last variable has only a diagonal entry */
    iRow[inz] = n-2;
    jCol[inz] = n-2;
    inz++;

    assert(inz == nele_hess);
  }
  else {
    /* return the values. This is a symmetric matrix, fill the upper
     * right triangle only */

    inz = 0;

    /* Diagonal entry for first variable */
    values[inz] = obj_factor*2.;
    inz++;

    for (i=1; i<n-1; i++) {
      values[inz] = obj_factor*2. + lambda[i-1]*2.*cos(x[i+1]);
      if (i>1) {
	values[inz] -= lambda[i-2]*(x[i-1]*x[i-1] + 1.5*x[i-1] -a[i-2])*cos(x[i]);
      }
      inz++;
      values[inz] = -lambda[i-1]*(2.*x[i]+1.5)*sin(x[i+1]);
      inz++;
    }

    values[inz] = obj_factor*2.;
    values[inz] -= lambda[n-3]*(x[n-2]*x[n-2] + 1.5*x[n-1] -a[n-3])*cos(x[n-1]);
    inz++;

    assert(inz == nele_hess);
  }

  return TRUE;
}
