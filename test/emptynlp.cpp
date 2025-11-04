// Copyright (C) 2021 COIN-OR Foundation
// All Rights Reserved.
// This code is published under the Eclipse Public License.

// get active asserts also if NDEBUG is defined
#ifdef NDEBUG
#undef NDEBUG
#endif

#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"
#include "IpTNLP.hpp"

#include <iostream>
#include <cassert>
#include <cmath>

using namespace Ipopt;

/** Empty NLP to test with Ipopt
 *
 * min sum_i x_i
 * s.t. x = 0           if not infeasbounds
 *      1 <= x <= 0     if infeasbounds
 *      sum_i x_i >= 0   if not infeascons and cons
 *      sum_i x_i >= 1   if infeascons
 */
class EmptyNLP: public TNLP
{
private:
   int nvars;
   bool cons;
   bool infeascons;
   bool infeasbounds;

public:
   /** constructor */
   EmptyNLP(
      int  nvars_         = 0,
      bool cons_          = true,
      bool infeascons_    = false,
      bool infeasbounds_  = false
   )
      : nvars(nvars_),
        cons(cons_),
        infeascons(infeascons_),
        infeasbounds(infeasbounds_)
   {
      assert(!infeasbounds || nvars > 0);
      assert(!infeascons || cons);
   }

   /** destructor */
   ~EmptyNLP() { }

   /** Method to return some info about the nlp */
   bool get_nlp_info(
      Index&          n,
      Index&          m,
      Index&          nnz_jac_g,
      Index&          nnz_h_lag,
      IndexStyleEnum& index_style
   )
   {
      n = nvars;
      m = cons ? 1 : 0;
      nnz_jac_g = cons ? n : 0;
      nnz_h_lag = 0;
      index_style = C_STYLE;

      return true;
   }

   /** Method to return the bounds for my problem */
   bool get_bounds_info(
      Index   n,
      Number* x_l,
      Number* x_u,
      Index   m,
      Number* g_l,
      Number* g_u
   )
   {
      for( Index i = 0; i < n; ++i )
      {
         x_l[i] = infeasbounds ? 1.0 : 0.0;
         x_u[i] = 0.0;
      }

      assert(m == (cons ? 1 : 0));
      if( cons )
      {
         g_l[0] = infeascons ? 1.0 : 0.0;
         g_u[0] = 1e300;
      }

      return true;
   }

   /** Method to return the starting point for the algorithm */
   bool get_starting_point(
      Index   n,
      bool    init_x,
      Number* x,
      bool    init_z,
      Number*,
      Number*,
      Index,
      bool    init_lambda,
      Number*
   )
   {
      if( init_x )
         for( Index i = 0; i < n; ++i )
         {
            x[i] = 10.0;
         }

      assert(!init_z);
      assert(!init_lambda);

      return true;
   }

   /** Method to return the objective value */
   bool eval_f(
      Index         n,
      const Number* x,
      bool,
      Number&       obj_value
   )
   {
      obj_value = 0.0;
      for( Index i = 0; i < n; ++i )
      {
         obj_value += x[i];
      }

      return true;
   }

   /** Method to return the gradient of the objective */
   bool eval_grad_f(
      Index         n,
      const Number*,
      bool,
      Number*       grad_f
   )
   {
      for( Index i = 0; i < n; ++i )
      {
         grad_f[i] = 1.0;
      }

      return true;
   }

   /** Method to return the constraint residuals */
   bool eval_g(
      Index         n,
      const Number* x,
      bool,
      Index         m,
      Number*       g
   )
   {
      if( !cons )
      {
         return true;
      }

      assert(m == 1);

      g[0] = 0.0;
      for( Index i = 0; i < n; ++i )
      {
         g[0] += x[i];
      }

      return true;
   }

   /** Method to return:
    *   1) The structure of the Jacobian (if "values" is NULL)
    *   2) The values of the Jacobian (if "values" is not NULL)
    */
   bool eval_jac_g(
      Index         n,
      const Number*,
      bool,
      Index,
      Index         nele_jac,
      Index*        iRow,
      Index*        jCol,
      Number*       values
   )
   {
      assert(nele_jac == (cons ? n : 0));

      if( !cons )
      {
         return true;
      }

      assert((iRow != NULL) == (jCol != NULL));
      assert((iRow != NULL) == (values == NULL));

      if( iRow != NULL )
         for( Index i = 0; i < n; ++i )
         {
            iRow[i] = 0;
            jCol[i] = i;
         }
      else
         for( Index i = 0; i < n; ++i )
         {
            values[i] = 1.0;
         }

      return true;
   }

   /** Method to return:
    *   1) The structure of the Hessian of the Lagrangian (if "values" is NULL)
    *   2) The values of the Hessian of the Lagrangian (if "values" is not NULL)
    */
   bool eval_h(
      Index,
      const Number*,
      bool,
      Number,
      Index,
      const Number*,
      bool,
      Index,
      Index*,
      Index*,
      Number*
   )
   {
      return true;
   }

   /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
   void finalize_solution(
      SolverReturn               status,
      Index                      n,
      const Number*              x,
      const Number*              z_L,
      const Number*              z_U,
      Index                      m,
      const Number*              g,
      const Number*              lambda,
      Number                     obj_value,
      const IpoptData*,
      IpoptCalculatedQuantities*
   )
   {
      std::cout << "Finalize called" << std::endl;
      std::cout << "x =";
      for( Index i = 0; i < n; ++i )
      {
         std::cout << ' ' << x[i];
      }
      std::cout << std::endl;
      std::cout << "z_L =";
      for( Index i = 0; i < n; ++i )
      {
         std::cout << ' ' << z_L[i];
      }
      std::cout << std::endl;
      std::cout << "z_U =";
      for( Index i = 0; i < n; ++i )
      {
         std::cout << ' ' << z_U[i];
      }
      std::cout << std::endl;
      std::cout << "lambda =";
      for( Index i = 0; i < m; ++i )
      {
         std::cout << ' ' << lambda[i];
      }
      std::cout << std::endl;

      if( status != (infeascons ? LOCAL_INFEASIBILITY : SUCCESS) )
      {
         std::cout << "*** Unexpected solver status " << status << " for" << (infeascons ? " infeasible" : "") << " NLP" << std::endl;
      }

      if( status == SUCCESS )
      {
         Number tol = 1e-6;
         assert(std::abs(obj_value) < tol);

         for( Index i = 0; i < n; ++i )
         {
            assert(std::abs(x[i]) < tol);
         }

         for( Index i = 0; i < m; ++i )
         {
            assert(std::abs(g[i]) < tol);
         }
      }
   }
};

bool runEmpty(
   int  nvars,
   bool cons,
   bool infeascons,
   bool infeasbounds
)
{
   std::cout << std::endl << "*** Solve for " << nvars << " variables, "
             << (cons ? 1 : 0) << ' '
             << (infeascons ? "infeasible" : "feasible") << " constraint, "
             << (infeasbounds ? "infeasible" : "feasible") << " bounds"
             << std::endl;

   // Create an instance of your nlp...
   SmartPtr<TNLP> nlp = new EmptyNLP(nvars, cons, infeascons, infeasbounds);

   // Create an instance of the IpoptApplication
   SmartPtr<IpoptApplication> app = new IpoptApplication();

   // Initialize the IpoptApplication and process the options
   ApplicationReturnStatus status;
   status = app->Initialize();
   if( status != Solve_Succeeded )
   {
      std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
      return false;
   }

   status = app->OptimizeTNLP(nlp);

   assert((status == Solve_Succeeded) == (!infeasbounds && !infeascons));
   assert((status == Infeasible_Problem_Detected) == infeascons);
   assert((status == Invalid_Problem_Definition) == infeasbounds);

   if( status > Not_Enough_Degrees_Of_Freedom )
   {
      assert(IsValid(app->Statistics()));

      // Retrieve some statistics about the solve
      Index iter_count = app->Statistics()->IterationCount();
      std::cout << std::endl << std::endl << "The problem solved in " << iter_count << " iterations!" << std::endl;

      Number final_obj = app->Statistics()->FinalObjective();
      std::cout << std::endl << std::endl << "The final value of the objective function is " << final_obj << '.'
                << std::endl;

      assert(app->Statistics()->IterationCount() == 0);
   }

   return true;
}

/** Mostly empty NLP to test reoptimization with Ipopt
 *
 * min sum_i x_i
 * s.t. 0 <= x
 *      x_i <= 0   for i <= nfixed
 *      sum_i x_i = rhs
 *      x_1 = 0    (as constraint)
 */
class ReOptNLP: public TNLP
{
private:
   Number rhs;
public:
   int nvars;
   int nfixed;

   /** constructor */
   ReOptNLP(
      Number rhs_ = 0.0
   )
      : rhs(rhs_),
        nvars(2),
        nfixed(0)
   { }

   /** Method to return some info about the nlp */
   bool get_nlp_info(
      Index&          n,
      Index&          m,
      Index&          nnz_jac_g,
      Index&          nnz_h_lag,
      IndexStyleEnum& index_style
   )
   {
      n = nvars;
      m = 2;
      nnz_jac_g = n + 1;
      nnz_h_lag = 0;
      index_style = C_STYLE;

      return true;
   }

   /** Method to return the bounds for my problem */
   bool get_bounds_info(
      Index   n,
      Number* x_l,
      Number* x_u,
      Index   m,
      Number* g_l,
      Number* g_u
   )
   {
      for( Index i = 0; i < n; ++i )
      {
         x_l[i] = 0.0;
         x_u[i] = i < nfixed ? 0.0 : 1.0;
      }

      assert(m == 2);
      g_l[0] = rhs;
      g_u[0] = rhs;

      g_l[1] = 0.0;
      g_u[1] = 0.0;

      return true;
   }

   /** Method to return the starting point for the algorithm */
   bool get_starting_point(
      Index   n,
      bool    init_x,
      Number* x,
      bool    init_z,
      Number*,
      Number*,
      Index,
      bool    init_lambda,
      Number*
   )
   {
      if( init_x )
         for( Index i = 0; i < n; ++i )
         {
            x[i] = 10.0;
         }

      assert(!init_z);
      assert(!init_lambda);

      return true;
   }

   /** Method to return the objective value */
   bool eval_f(
      Index         n,
      const Number* x,
      bool,
      Number&       obj_value
   )
   {
      obj_value = 0.0;
      for( Index i = 0; i < n; ++i )
      {
         obj_value += x[i];
      }

      return true;
   }

   /** Method to return the gradient of the objective */
   bool eval_grad_f(
      Index         n,
      const Number*,
      bool,
      Number*       grad_f
   )
   {
      for( Index i = 0; i < n; ++i )
      {
         grad_f[i] = 1.0;
      }

      return true;
   }

   /** Method to return the constraint residuals */
   bool eval_g(
      Index         n,
      const Number* x,
      bool,
      Index         m,
      Number*       g
   )
   {
      assert(m == 2);

      g[0] = 0.0;
      for( Index i = 0; i < n; ++i )
      {
         g[0] += x[i];
      }

      g[1] = x[0];

      return true;
   }

   /** Method to return:
    *   1) The structure of the Jacobian (if "values" is NULL)
    *   2) The values of the Jacobian (if "values" is not NULL)
    */
   bool eval_jac_g(
      Index         n,
      const Number*,
      bool,
      Index,
      Index         nele_jac,
      Index*        iRow,
      Index*        jCol,
      Number*       values
   )
   {
      assert((iRow != NULL) == (jCol != NULL));
      assert((iRow != NULL) == (values == NULL));
      assert(nele_jac == n + 1);

      if( iRow != NULL )
      {
         for( Index i = 0; i < n; ++i )
         {
            iRow[i] = 0;
            jCol[i] = i;
         }
         iRow[n] = 1;
         jCol[n] = 0;
      }
      else
         for( Index i = 0; i < n + 1; ++i )
         {
            values[i] = 1.0;
         }

      return true;
   }

   /** Method to return:
    *   1) The structure of the Hessian of the Lagrangian (if "values" is NULL)
    *   2) The values of the Hessian of the Lagrangian (if "values" is not NULL)
    */
   bool eval_h(
      Index,
      const Number*,
      bool,
      Number,
      Index,
      const Number*,
      bool,
      Index,
      Index*,
      Index*,
      Number*
   )
   {
      return true;
   }

   /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
   void finalize_solution(
      SolverReturn               status,
      Index                      n,
      const Number*              x,
      const Number*              z_L,
      const Number*              z_U,
      Index                      m,
      const Number*              g,
      const Number*              lambda,
      Number                     obj_value,
      const IpoptData*,
      IpoptCalculatedQuantities*
   )
   {
      std::cout << "Finalize called" << std::endl;
      std::cout << "x =";
      for( Index i = 0; i < n; ++i )
      {
         std::cout << ' ' << x[i];
      }
      std::cout << std::endl;
      std::cout << "z_L =";
      for( Index i = 0; i < n; ++i )
      {
         std::cout << ' ' << z_L[i];
      }
      std::cout << std::endl;
      std::cout << "z_U =";
      for( Index i = 0; i < n; ++i )
      {
         std::cout << ' ' << z_U[i];
      }
      std::cout << std::endl;
      std::cout << "lambda =";
      for( Index i = 0; i < m; ++i )
      {
         std::cout << ' ' << lambda[i];
      }
      std::cout << std::endl;

      if( status == SUCCESS )
      {
         Number tol = 1e-5;
         assert(std::abs(obj_value - rhs) < tol);

         if( rhs == 0.0 )
         {
            for( Index i = 0; i < n; ++i )
            {
               assert(std::abs(x[i]) < tol);
            }
         }

         assert(std::abs(g[0] - rhs) < tol);
         assert(std::abs(g[1]) < tol);
      }
   }
};

bool runReOpt(
   int nvars1,
   int nvars2,
   int nfixed,
   Number rhs
)
{
   std::cout << std::endl << "*** Solve with " << nvars1 << " variables and rhs=" << rhs << "." << std::endl;

   SmartPtr<ReOptNLP> nlp = new ReOptNLP(rhs);
   nlp->nvars = nvars1;
   SmartPtr<IpoptApplication> app = new IpoptApplication();

   ApplicationReturnStatus status;
   status = app->Initialize();
   if( status != Solve_Succeeded )
   {
      std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
      return false;
   }

   status = app->OptimizeTNLP(GetRawPtr(nlp));

   if( nvars1 < 2 )
   {
      assert(status == Not_Enough_Degrees_Of_Freedom);
   }
   else
   {
      assert((status == Solve_Succeeded) == (rhs <= 1));
      assert((status == Infeasible_Problem_Detected) == (rhs > 1));
   }

   std::cout << std::endl << "*** Resolve with " << nvars2 << " variables, " << nfixed << " variables fixed." << std::endl;

   nlp->nvars = nvars2;
   nlp->nfixed = nfixed;
   status = app->ReOptimizeTNLP(GetRawPtr(nlp));

   if( nvars2 < 2 && nfixed < nvars2 )
   {
      assert(status == Not_Enough_Degrees_Of_Freedom);
   }
   else
   {
      assert((status == Solve_Succeeded) == (rhs <= 1 && (rhs == 0.0 || nfixed < nvars2)));
      assert((status == Infeasible_Problem_Detected) == (rhs > 1 || (rhs != 0.0 && nfixed == nvars2)));
   }

   return true;
}

int main(
   int,
   char**
)
{
   if( !runEmpty(0, false, false, false) )
   {
      return EXIT_FAILURE;
   }

   if( !runEmpty(0, true, false, false) )
   {
      return EXIT_FAILURE;
   }

   if( !runEmpty(5, true, false, false) )
   {
      return EXIT_FAILURE;
   }

   if( !runEmpty(0, true, true, false) )
   {
      return EXIT_FAILURE;
   }

   if( !runEmpty(5, true, true, false) )
   {
      return EXIT_FAILURE;
   }

   if( !runEmpty(5, true, false, true) )
   {
      return EXIT_FAILURE;
   }

   // 2 variables, fix them in resolve, 2 constraints
   if( !runReOpt(2, 2, 2, 0.0) )
   {
      return EXIT_FAILURE;
   }

   if( !runReOpt(2, 2, 2, 1.0) )
   {
      return EXIT_FAILURE;
   }

   // 1 variable, do not fix in resolve, 2 constraints
   // should give a too-few-degree-of-freedom in both solves
   if( !runReOpt(1, 1, 0, 0.0) )
   {
      return EXIT_FAILURE;
   }

   // 2 variables, then 1 variable, do not fix in resolve
   // should give a too-few-degree-of-freedom for second solve  (this case had a bug)
   if( !runReOpt(2, 1, 0, 0.0) )
   {
      return EXIT_FAILURE;
   }

   // 2 variables, then 1 variable, but fixed in resolve
   // should not give a too-few-degree-of-freedom, but feasible
   if( !runReOpt(2, 1, 1, 0.0) )
   {
      return EXIT_FAILURE;
   }

   // 2 variables, then 1 variable, but fixed in resolve
   // should not give a too-few-degree-of-freedom, but infeasible
   if( !runReOpt(2, 1, 1, 1.0) )
   {
      return EXIT_FAILURE;
   }

   std::cout << std::endl << "*** All tests passed" << std::endl;

   return EXIT_SUCCESS;
}
