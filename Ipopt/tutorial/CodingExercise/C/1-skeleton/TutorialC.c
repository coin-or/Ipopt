/* Copyright (C) 2009 International Business Machines.
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * $Id$
 *
 * Author:  Andreas Waechter               IBM    2009-04-02
 */

/*
// This file is part of the Ipopt tutorial.  It is the skeleton for
// the C implemention of the coding exercise problem (in AMPL
// formulation):
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
  size = 300;

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
  SET BOUNDS

  /* set the number of constraints and allocate space for the bounds */
  m=size-2;
  g_L = (Number*)malloc(sizeof(Number)*m);
  g_U = (Number*)malloc(sizeof(Number)*m);
  /* set the values of the constraint bounds */
  SET BOUNDS
  
  /* Number of nonzeros in the Jacobian of the constraints
     each constraint has three nonzeros */
  nele_jac = FILLME

  /* Number of nonzeros in the Hessian of the Lagrangian (lower or
     upper triangual part only)
     We have the full diagonal, and the first off-diagonal except for
     the first and last variable */
  nele_hess = FILLME

  /* indexing style for matrices */
  index_style = FILLME; /* C-style; start counting of rows and column
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
  SET INITIAL POINT

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
  double* a = PD->a;

  return TRUE;
}

Bool eval_grad_f(Index n, Number* x, Bool new_x,
                 Number* grad_f, UserDataPtr user_data)
{
  int i;
  ProblemData PD = (ProblemData)user_data;
  double* a = PD->a;

  return TRUE;
}

Bool eval_g(Index n, Number* x, Bool new_x,
            Index m, Number* g, UserDataPtr user_data)
{
  int j;
  ProblemData PD = (ProblemData)user_data;
  double* a = PD->a;

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


  return TRUE;
}
