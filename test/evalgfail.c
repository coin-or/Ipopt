/* Copyright (C) 2026 Steven R. Hall
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * A constraint callback that reports failure once, with check_derivatives_for_naninf on.
 *
 * Hock-Schittkowski problem 71 through the C interface, as in examples/hs071_c, except that
 * eval_g returns false on its third call, a trial point of the first line search. Ipopt must
 * treat that as an evaluation error, backtrack, and still find the optimal solution. It used
 * to take the norm of the constraint vector, which the failed evaluation never filled, before
 * deciding whether to report an invalid number, and crash on the unallocated values.
 */

#include "IpStdCInterface.h"
#include <stdio.h>
#include <stdlib.h>

static int g_calls = 0;

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
   if( ++g_calls == 3 )
   {
      printf("eval_g: reporting failure on call %d\n", g_calls);
      return false;
   }
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

int main(void)
{
   ipnumber x_L[4] = { 1.0, 1.0, 1.0, 1.0 };
   ipnumber x_U[4] = { 5.0, 5.0, 5.0, 5.0 };
   ipnumber g_L[2] = { 25.0, 40.0 };
   ipnumber g_U[2] = { 2e19, 40.0 };
   ipnumber x[4] = { 1.0, 5.0, 5.0, 1.0 };
   ipnumber obj;
   IpoptProblem nlp;
   enum ApplicationReturnStatus status;

   nlp = CreateIpoptProblem(4, x_L, x_U, 2, g_L, g_U, 8, 0, 0,
                            &eval_f, &eval_g, &eval_grad_f, &eval_jac_g, &eval_h);
   AddIpoptStrOption(nlp, "hessian_approximation", "limited-memory");
   AddIpoptStrOption(nlp, "check_derivatives_for_naninf", "yes");

   status = IpoptSolve(nlp, x, NULL, &obj, NULL, NULL, NULL, NULL);
   FreeIpoptProblem(nlp);

   if( status != Solve_Succeeded )
   {
      printf("IpoptSolve returned %d\n", (int) status);
      return EXIT_FAILURE;
   }
   return EXIT_SUCCESS;
}
