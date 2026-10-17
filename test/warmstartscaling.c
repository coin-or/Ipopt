/* Copyright (C) 2026 Steven R. Hall
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * Warm-start bound multipliers under user scaling of the variables.
 *
 * Hock-Schittkowski problem 71 through the C interface, as in examples/hs071_c, with
 * nlp_scaling_method = user-scaling and a scaling factor other than 1 on x[0], whose lower bound
 * is active at the solution. The problem is solved, then solved again from that solution and
 * its multipliers with warm_start_init_point = yes and max_iter = 0, so the multipliers
 * returned are those Ipopt started from. They must equal the ones it was given. The initial
 * bound multipliers used to be scaled by s_x where the inverse of the conversion of the final
 * ones needs 1 / s_x, so they came back multiplied by s_x^2.
 */

#include "IpStdCInterface.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static bool eval_f(
   ipindex     n,
   ipnumber*   x,
   bool        new_x,
   ipnumber*   obj_value,
   UserDataPtr user_data
)
{
   (void) n;
   (void) new_x;
   (void) user_data;
   *obj_value = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
   return true;
}

static bool eval_grad_f(
   ipindex     n,
   ipnumber*   x,
   bool        new_x,
   ipnumber*   grad_f,
   UserDataPtr user_data
)
{
   (void) n;
   (void) new_x;
   (void) user_data;
   grad_f[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
   grad_f[1] = x[0] * x[3];
   grad_f[2] = x[0] * x[3] + 1;
   grad_f[3] = x[0] * (x[0] + x[1] + x[2]);
   return true;
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
   (void) n;
   (void) new_x;
   (void) m;
   (void) user_data;
   g[0] = x[0] * x[1] * x[2] * x[3];
   g[1] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3];
   return true;
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
   (void) n;
   (void) new_x;
   (void) m;
   (void) nele_jac;
   (void) user_data;
   if( values == NULL )
   {
      /* dense 2 by 4 */
      int k = 0;
      for( int i = 0; i < 2; i++ )
      {
         for( int j = 0; j < 4; j++ )
         {
            iRow[k] = i;
            jCol[k] = j;
            k++;
         }
      }
      return true;
   }
   values[0] = x[1] * x[2] * x[3];
   values[1] = x[0] * x[2] * x[3];
   values[2] = x[0] * x[1] * x[3];
   values[3] = x[0] * x[1] * x[2];
   values[4] = 2 * x[0];
   values[5] = 2 * x[1];
   values[6] = 2 * x[2];
   values[7] = 2 * x[3];
   return true;
}

/* never called: the Hessian is approximated, but the interface requires a callback */
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
   (void) n;
   (void) x;
   (void) new_x;
   (void) obj_factor;
   (void) m;
   (void) lambda;
   (void) new_lambda;
   (void) nele_hess;
   (void) iRow;
   (void) jCol;
   (void) values;
   (void) user_data;
   return true;
}

static IpoptProblem create_problem(void)
{
   ipnumber x_L[4] = { 1.0, 1.0, 1.0, 1.0 };
   ipnumber x_U[4] = { 5.0, 5.0, 5.0, 5.0 };
   ipnumber g_L[2] = { 25.0, 40.0 };
   ipnumber g_U[2] = { 2e19, 40.0 };
   ipnumber x_scaling[4] = { 2.0, 1.0, 1.0, 1.0 };
   ipnumber g_scaling[2] = { 0.5, 0.5 };
   IpoptProblem nlp;

   nlp = CreateIpoptProblem(4, x_L, x_U, 2, g_L, g_U, 8, 0, 0,
                            &eval_f, &eval_g, &eval_grad_f, &eval_jac_g, &eval_h);
   AddIpoptStrOption(nlp, "hessian_approximation", "limited-memory");
   AddIpoptStrOption(nlp, "nlp_scaling_method", "user-scaling");
   SetIpoptProblemScaling(nlp, 3.0, x_scaling, g_scaling);
   AddIpoptIntOption(nlp, "print_level", 0);
   AddIpoptStrOption(nlp, "sb", "yes");
   return nlp;
}

static bool agrees(
   const char* name,
   ipindex     i,
   ipnumber    given,
   ipnumber    returned
)
{
   if( fabs(returned - given) > 1e-6 * (1.0 + fabs(given)) )
   {
      printf("%s[%d]: given %g, returned %g\n", name, (int) i, given, returned);
      return false;
   }
   return true;
}

int main(void)
{
   ipnumber x[4] = { 1.0, 5.0, 5.0, 1.0 };
   ipnumber obj;
   ipnumber mult_g[2];
   ipnumber mult_x_L[4];
   ipnumber mult_x_U[4];
   ipnumber given_g[2];
   ipnumber given_L[4];
   ipnumber given_U[4];
   IpoptProblem nlp;
   enum ApplicationReturnStatus status;
   bool ok = true;
   ipindex i;

   nlp = create_problem();
   status = IpoptSolve(nlp, x, NULL, &obj, mult_g, mult_x_L, mult_x_U, NULL);
   FreeIpoptProblem(nlp);
   if( status != Solve_Succeeded )
   {
      printf("IpoptSolve returned %d\n", (int) status);
      return EXIT_FAILURE;
   }

   for( i = 0; i < 2; ++i )
   {
      given_g[i] = mult_g[i];
   }
   for( i = 0; i < 4; ++i )
   {
      given_L[i] = mult_x_L[i];
      given_U[i] = mult_x_U[i];
   }

   /* start from the solution and its multipliers, without moving them */
   nlp = create_problem();
   AddIpoptStrOption(nlp, "warm_start_init_point", "yes");
   AddIpoptNumOption(nlp, "warm_start_bound_push", 1e-9);
   AddIpoptNumOption(nlp, "warm_start_bound_frac", 1e-9);
   AddIpoptNumOption(nlp, "warm_start_slack_bound_push", 1e-9);
   AddIpoptNumOption(nlp, "warm_start_slack_bound_frac", 1e-9);
   AddIpoptNumOption(nlp, "warm_start_mult_bound_push", 1e-9);
   AddIpoptIntOption(nlp, "max_iter", 0);
   IpoptSolve(nlp, x, NULL, &obj, mult_g, mult_x_L, mult_x_U, NULL);
   FreeIpoptProblem(nlp);

   for( i = 0; i < 2; ++i )
   {
      ok = agrees("mult_g", i, given_g[i], mult_g[i]) && ok;
   }
   for( i = 0; i < 4; ++i )
   {
      ok = agrees("mult_x_L", i, given_L[i], mult_x_L[i]) && ok;
      ok = agrees("mult_x_U", i, given_U[i], mult_x_U[i]) && ok;
   }
   return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
