/* Copyright (C) 2005, 2011 International Business Machines and others.
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-17
 */

#include "IpStdCInterface.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

/** This is an example how user_data can be used. */
struct MyUserData
{
   ipnumber g_offset[2]; /**< This is an offset for the constraints.  */
   IpoptProblem nlp;   /**< The problem to be solved. Required in intermediate_cb. */
};

/* Callback Implementations */
static
bool eval_f(
   ipindex     n,
   ipnumber*   x,
   bool        new_x,
   ipnumber*   obj_value,
   UserDataPtr user_data
)
{
   assert(n == 4);
   (void) n;

   (void) new_x;
   (void) user_data;

   *obj_value = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];

   return true;
}

static
bool eval_grad_f(
   ipindex     n,
   ipnumber*   x,
   bool        new_x,
   ipnumber*   grad_f,
   UserDataPtr user_data
)
{
   assert(n == 4);
   (void) n;

   (void) new_x;
   (void) user_data;

   grad_f[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
   grad_f[1] = x[0] * x[3];
   grad_f[2] = x[0] * x[3] + 1;
   grad_f[3] = x[0] * (x[0] + x[1] + x[2]);

   return true;
}

static
bool eval_g(
   ipindex     n,
   ipnumber*   x,
   bool        new_x,
   ipindex     m,
   ipnumber*   g,
   UserDataPtr user_data
)
{
   struct MyUserData* my_data = user_data;

   assert(n == 4);
   (void) n;
   assert(m == 2);
   (void) m;

   (void) new_x;

   g[0] = x[0] * x[1] * x[2] * x[3] + my_data->g_offset[0];
   g[1] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] + my_data->g_offset[1];

   return true;
}

static
bool eval_jac_g(
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
   (void) n;
   (void) new_x;
   (void) m;
   (void) nele_jac;
   (void) user_data;

   if( values == NULL )
   {
      /* return the structure of the jacobian */

      /* this particular jacobian is dense */
      iRow[0] = 0;
      jCol[0] = 0;
      iRow[1] = 0;
      jCol[1] = 1;
      iRow[2] = 0;
      jCol[2] = 2;
      iRow[3] = 0;
      jCol[3] = 3;
      iRow[4] = 1;
      jCol[4] = 0;
      iRow[5] = 1;
      jCol[5] = 1;
      iRow[6] = 1;
      jCol[6] = 2;
      iRow[7] = 1;
      jCol[7] = 3;
   }
   else
   {
      /* return the values of the jacobian of the constraints */

      values[0] = x[1] * x[2] * x[3]; /* 0,0 */
      values[1] = x[0] * x[2] * x[3]; /* 0,1 */
      values[2] = x[0] * x[1] * x[3]; /* 0,2 */
      values[3] = x[0] * x[1] * x[2]; /* 0,3 */

      values[4] = 2 * x[0]; /* 1,0 */
      values[5] = 2 * x[1]; /* 1,1 */
      values[6] = 2 * x[2]; /* 1,2 */
      values[7] = 2 * x[3]; /* 1,3 */
   }

   return true;
}

static
bool eval_h(
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
   (void) n;
   (void) new_x;
   (void) m;
   (void) new_lambda;
   (void) user_data;

   if( values == NULL )
   {
      ipindex idx; /* nonzero element counter */
      ipindex row; /* row counter for loop */
      ipindex col; /* col counter for loop */

      /* return the structure. This is a symmetric matrix, fill the lower left
       * triangle only. */

      /* the hessian for this problem is actually dense */
      idx = 0;
      for( row = 0; row < 4; row++ )
      {
         for( col = 0; col <= row; col++ )
         {
            iRow[idx] = row;
            jCol[idx] = col;
            idx++;
         }
      }

      assert(idx == nele_hess);
      (void) nele_hess;
   }
   else
   {
      /* return the values. This is a symmetric matrix, fill the lower left
       * triangle only */

      /* fill the objective portion */
      values[0] = obj_factor * (2 * x[3]); /* 0,0 */

      values[1] = obj_factor * (x[3]); /* 1,0 */
      values[2] = 0; /* 1,1 */

      values[3] = obj_factor * (x[3]); /* 2,0 */
      values[4] = 0; /* 2,1 */
      values[5] = 0; /* 2,2 */

      values[6] = obj_factor * (2 * x[0] + x[1] + x[2]); /* 3,0 */
      values[7] = obj_factor * (x[0]); /* 3,1 */
      values[8] = obj_factor * (x[0]); /* 3,2 */
      values[9] = 0; /* 3,3 */

      /* add the portion for the first constraint */
      values[1] += lambda[0] * (x[2] * x[3]); /* 1,0 */

      values[3] += lambda[0] * (x[1] * x[3]); /* 2,0 */
      values[4] += lambda[0] * (x[0] * x[3]); /* 2,1 */

      values[6] += lambda[0] * (x[1] * x[2]); /* 3,0 */
      values[7] += lambda[0] * (x[0] * x[2]); /* 3,1 */
      values[8] += lambda[0] * (x[0] * x[1]); /* 3,2 */

      /* add the portion for the second constraint */
      values[0] += lambda[1] * 2; /* 0,0 */

      values[2] += lambda[1] * 2; /* 1,1 */

      values[5] += lambda[1] * 2; /* 2,2 */

      values[9] += lambda[1] * 2; /* 3,3 */
   }

   return true;
}

static IPOPT_UNUSED
bool intermediate_cb(
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
   bool have_iter;
   bool have_viol;

   ipnumber x[4];
   ipnumber x_L_viol[4];
   ipnumber x_U_viol[4];
   ipnumber z_L[4];
   ipnumber z_U[4];
   ipnumber compl_x_L[4];
   ipnumber compl_x_U[4];
   ipnumber grad_lag_x[4];

   ipnumber g[2];
   ipnumber lambda[2];
   ipnumber constraint_violation[2];
   ipnumber compl_g[2];

   IpoptProblem nlp = ((struct MyUserData*)user_data)->nlp;

   have_iter = GetIpoptCurrentIterate(nlp, false, 4, x, z_L, z_U, 2, g, lambda);
   have_viol = GetIpoptCurrentViolations(nlp, false, 4, x_L_viol, x_U_viol, compl_x_L, compl_x_U, grad_lag_x, 2, constraint_violation, compl_g);

   printf("Current iterate at iteration %d:\n", (int)iter_count);
   printf("  %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n", "x", "z_L", "z_U", "bound_viol", "compl_x_L", "compl_x_U", "grad_lag_x");
   for( int i = 0; i < 4; ++i )
   {
      if( have_iter )
      {
         printf("  %-12g %-12g %-12g", x[i], z_L[i], z_U[i]);
      }
      else
      {
         printf("  %-12s %-12s %-12s", "n/a", "n/a", "n/a");
      }
      if( have_viol )
      {
         printf(" %-12g %-12g %-12g %-12g\n", x_L_viol[i] > x_U_viol[i] ? x_L_viol[i] : x_U_viol[i], compl_x_L[i], compl_x_U[i], grad_lag_x[i]);
      }
      else
      {
         printf(" %-12s %-12s %-12s %-12s\n", "n/a", "n/a", "n/a", "n/a");
      }
   }

   printf("  %-12s %-12s %-12s %-12s\n", "g(x)", "lambda", "constr_viol", "compl_g");
   for( int i = 0; i < 2; ++i )
   {
      if( have_iter )
      {
         printf("  %-12g %-12g", g[i], lambda[i]);
      }
      else
      {
         printf("  %-12s %-12s", "n/a", "n/a");
      }
      if( have_viol )
      {
         printf(" %-12g %-12g\n", constraint_violation[i], compl_g[i]);
      }
      else
      {
         printf(" %-12s %-12s\n", "n/a", "n/a");
      }
   }

   if( inf_pr < 1e-4 )
   {
      return false;
   }

   (void) alg_mod;
   (void) obj_value;
   (void) inf_du;
   (void) mu;
   (void) d_norm;
   (void) regularization_size;
   (void) alpha_du;
   (void) alpha_pr;
   (void) ls_trials;

   return true;
}

/** Main Program */
/* [MAIN] */
int main()
{
   ipindex n = -1;                      /* number of variables */
   ipindex m = -1;                      /* number of constraints */
   ipindex nele_jac;                    /* number of nonzeros in the Jacobian of the constraints */
   ipindex nele_hess;                   /* number of nonzeros in the Hessian of the Lagrangian (lower or upper triangular part only) */
   ipindex index_style;                 /* indexing style for matrices */
   ipnumber* x_L = NULL;                /* lower bounds on x */
   ipnumber* x_U = NULL;                /* upper bounds on x */
   ipnumber* g_L = NULL;                /* lower bounds on g */
   ipnumber* g_U = NULL;                /* upper bounds on g */
   IpoptProblem nlp = NULL;             /* IpoptProblem */
   enum ApplicationReturnStatus status; /* Solve return code */
   ipnumber* x = NULL;                  /* starting point and solution vector */
   ipnumber* mult_g = NULL;             /* constraint multipliers at the solution */
   ipnumber* mult_x_L = NULL;           /* lower bound multipliers at the solution */
   ipnumber* mult_x_U = NULL;           /* upper bound multipliers at the solution */
   ipnumber obj;                        /* objective value */
   struct MyUserData user_data;         /* our user data for the function evaluations */
   ipindex i;                           /* generic counter */

   /* set the number of variables and allocate space for the bounds */
   n = 4;
   x_L = (ipnumber*) malloc(sizeof(ipnumber) * n);
   x_U = (ipnumber*) malloc(sizeof(ipnumber) * n);
   /* set the values for the variable bounds */
   for( i = 0; i < n; i++ )
   {
      x_L[i] = 1.0;
      x_U[i] = 5.0;
   }

   /* set the number of constraints and allocate space for the bounds */
   m = 2;
   g_L = (ipnumber*) malloc(sizeof(ipnumber) * m);
   g_U = (ipnumber*) malloc(sizeof(ipnumber) * m);
   /* set the values of the constraint bounds */
   g_L[0] = 25;
   g_U[0] = 2e19;
   g_L[1] = 40;
   g_U[1] = 40;

   /* set the number of nonzeros in the Jacobian and Hessian */
   nele_jac = 8;
   nele_hess = 10;

   /* set the indexing style to C-style (start counting of rows and column indices at 0) */
   index_style = 0;

   /* create the IpoptProblem */
   nlp = CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess, index_style,
                            &eval_f, &eval_g, &eval_grad_f,
                            &eval_jac_g, &eval_h);

   /* We can free the memory now - the values for the bounds have been
    * copied internally in CreateIpoptProblem
    */
   free(x_L);
   free(x_U);
   free(g_L);
   free(g_U);

   /* Set some options.  Note the following ones are only examples,
    * they might not be suitable for your problem.
    */
   AddIpoptNumOption(nlp, "tol", 3.82e-6);
   AddIpoptStrOption(nlp, "mu_strategy", "adaptive");
   AddIpoptStrOption(nlp, "output_file", "ipopt.out");

   /* allocate space for the initial point and set the values */
   x = (ipnumber*) malloc(sizeof(ipnumber) * n);
   x[0] = 1.0;
   x[1] = 5.0;
   x[2] = 5.0;
   x[3] = 1.0;

   /* allocate space to store the bound multipliers at the solution */
   mult_g = (ipnumber*) malloc(sizeof(ipnumber) * m);
   mult_x_L = (ipnumber*) malloc(sizeof(ipnumber) * n);
   mult_x_U = (ipnumber*) malloc(sizeof(ipnumber) * n);

   /* Initialize the user data */
   user_data.g_offset[0] = 0.;
   user_data.g_offset[1] = 0.;
   user_data.nlp = nlp;

   /* Set the callback method for intermediate user-control.
    * This is not required, just gives you some intermediate control in
    * case you need it.
    */
   /* SetIntermediateCallback(nlp, intermediate_cb); */

   /* solve the problem */
   status = IpoptSolve(nlp, x, NULL, &obj, mult_g, mult_x_L, mult_x_U, &user_data);

   if( status == Solve_Succeeded )
   {
      printf("\n\nSolution of the primal variables, x\n");
      for( i = 0; i < n; i++ )
      {
         printf("x[%d] = %e\n", (int)i, x[i]);
      }

      printf("\n\nSolution of the constraint multipliers, lambda\n");
      for( i = 0; i < m; i++ )
      {
         printf("lambda[%d] = %e\n", (int)i, mult_g[i]);
      }
      printf("\n\nSolution of the bound multipliers, z_L and z_U\n");
      for( i = 0; i < n; i++ )
      {
         printf("z_L[%d] = %e\n", (int)i, mult_x_L[i]);
      }
      for( i = 0; i < n; i++ )
      {
         printf("z_U[%d] = %e\n", (int)i, mult_x_U[i]);
      }

      printf("\n\nObjective value\nf(x*) = %e\n", obj);
   }
   else
   {
      printf("\n\nERROR OCCURRED DURING IPOPT OPTIMIZATION.\n");
   }

   /* Now we are going to solve this problem again, but with slightly
    * modified constraints.  We change the constraint offset of the
    * first constraint a bit, and resolve the problem using the warm
    * start option.
    */
   user_data.g_offset[0] = 0.2;

   if( status == Solve_Succeeded )
   {
      /* Now resolve with a warmstart. */
      AddIpoptStrOption(nlp, "warm_start_init_point", "yes");
      /* The following option reduce the automatic modification of the
       * starting point done my Ipopt.
       */
      AddIpoptNumOption(nlp, "bound_push", 1e-5);
      AddIpoptNumOption(nlp, "bound_frac", 1e-5);
      status = IpoptSolve(nlp, x, NULL, &obj, mult_g, mult_x_L, mult_x_U, &user_data);

      if( status == Solve_Succeeded )
      {
         printf("\n\nSolution of the primal variables, x\n");
         for( i = 0; i < n; i++ )
         {
            printf("x[%d] = %e\n", (int)i, x[i]);
         }

         printf("\n\nSolution of the constraint multipliers, lambda\n");
         for( i = 0; i < m; i++ )
         {
            printf("lambda[%d] = %e\n", (int)i, mult_g[i]);
         }
         printf("\n\nSolution of the bound multipliers, z_L and z_U\n");
         for( i = 0; i < n; i++ )
         {
            printf("z_L[%d] = %e\n", (int)i, mult_x_L[i]);
         }
         for( i = 0; i < n; i++ )
         {
            printf("z_U[%d] = %e\n", (int)i, mult_x_U[i]);
         }

         printf("\n\nObjective value\nf(x*) = %e\n", obj);
      }
      else
      {
         printf("\n\nERROR OCCURRED DURING IPOPT OPTIMIZATION WITH WARM START.\n");
      }
   }

   /* free allocated memory */
   FreeIpoptProblem(nlp);
   free(x);
   free(mult_g);
   free(mult_x_L);
   free(mult_x_U);

   return (status == Solve_Succeeded) ? EXIT_SUCCESS : EXIT_FAILURE;
}
/* [MAIN] */
