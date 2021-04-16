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
 *      sum_i x_i >= 0   if not infeascons
 *      sum_i x_i >= 1   if infeascons
 */
class EmptyNLP: public TNLP
{
private:
   int nvars;
   bool infeascons;
   bool infeasbounds;

public:
   /** constructor */
   EmptyNLP(
      int  nvars_         = 0,
      bool infeascons_    = false,
      bool infeasbounds_  = false
   )
   : nvars(nvars_),
     infeascons(infeascons_),
     infeasbounds(infeasbounds_)
   {
      assert(!infeasbounds || nvars > 0);
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
      m = 1;
      nnz_jac_g = n;
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

      assert(m == 1);
      g_l[0] = infeascons ? 1.0 : 0.0;
      g_u[0] = 1e300;

      return true;
   }

   /** Method to return the starting point for the algorithm */
   bool get_starting_point(
      Index   n,
      bool    init_x,
      Number* x,
      bool    init_z,
      Number* ,
      Number* ,
      Index   ,
      bool    init_lambda,
      Number*
   )
   {
      if( init_x )
         for( Index i = 0; i < n; ++i )
            x[i] = 10.0;

      assert(!init_z);
      assert(!init_lambda);

      return true;
   }

   /** Method to return the objective value */
   bool eval_f(
      Index         n,
      const Number* x,
      bool          ,
      Number&       obj_value
   )
   {
      obj_value = 0.0;
      for( int i = 0; i < n; ++i )
         obj_value += x[i];

      return true;
   }

   /** Method to return the gradient of the objective */
   bool eval_grad_f(
      Index         n,
      const Number* ,
      bool          ,
      Number*       grad_f
   )
   {
      for( int i = 0; i < n; ++i )
         grad_f[i] = 1.0;

      return true;
   }

   /** Method to return the constraint residuals */
   bool eval_g(
      Index         n,
      const Number* x,
      bool          ,
      Index         m,
      Number*       g
   )
   {
      assert(m == 1);

      g[0] = 0.0;
      for( int i = 0; i < n; ++i )
         g[0] += x[i];

      return true;
   }

   /** Method to return:
    *   1) The structure of the Jacobian (if "values" is NULL)
    *   2) The values of the Jacobian (if "values" is not NULL)
    */
   bool eval_jac_g(
      Index         n,
      const Number* ,
      bool          ,
      Index         ,
      Index         nele_jac,
      Index*        iRow,
      Index*        jCol,
      Number*       values
   )
   {
      assert((iRow != NULL) == (jCol != NULL));
      assert((iRow != NULL) == (values == NULL));
      assert(nele_jac == n);

      if( iRow != NULL )
         for( int i = 0; i < n; ++i )
         {
            iRow[i] = 0;
            jCol[i] = i;
         }
      else
         for( int i = 0; i < n; ++i )
            values[i] = 1.0;

      return true;
   }

   /** Method to return:
    *   1) The structure of the Hessian of the Lagrangian (if "values" is NULL)
    *   2) The values of the Hessian of the Lagrangian (if "values" is not NULL)
    */
   bool eval_h(
      Index         ,
      const Number* ,
      bool          ,
      Number        ,
      Index         ,
      const Number* ,
      bool          ,
      Index         ,
      Index*        ,
      Index*        ,
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
      const IpoptData*           ,
      IpoptCalculatedQuantities*
   )
   {
      std::cout << "Finalize called" << std::endl;
      std::cout << "x =";
      for( int i = 0; i < n; ++i )
         std::cout << ' ' << x[i];
      std::cout << std::endl;
      std::cout << "z_L =";
      for( int i = 0; i < n; ++i )
         std::cout << ' ' << z_L[i];
      std::cout << std::endl;
      std::cout << "z_U =";
      for( int i = 0; i < n; ++i )
         std::cout << ' ' << z_U[i];
      std::cout << std::endl;
      std::cout << "lambda =";
      for( int i = 0; i < m; ++i )
         std::cout << ' ' << lambda[i];
      std::cout << std::endl;

      if( status != (infeascons ? LOCAL_INFEASIBILITY : SUCCESS) )
      {
         std::cout << "*** Unexpected solver status " << status << " for" << (infeascons ? " infeasible" : "") << " NLP" << std::endl;
      }

      if( status == SUCCESS )
      {
         Number tol = 1e-6;
         assert(fabs(obj_value) < tol);

         for( int i = 0; i < n; ++i )
            assert(fabs(x[i]) < tol);

         for( int i = 0; i < m; ++i )
            assert(fabs(g[i]) < tol);
      }
   }
};

bool run(
   int  nvars,
   bool infeascons,
   bool infeasbounds
   )
{
   std::cout << std::endl << "*** Solve for " << nvars << " variables, "
      << (infeascons ? "infeasible" : "feasible") << " constraint, "
      << (infeasbounds ? "infeasible" : "feasible") << " bounds"
      << std::endl;

   // Create an instance of your nlp...
   SmartPtr<TNLP> nlp = new EmptyNLP(nvars, infeascons, infeasbounds);

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

int main(
   int,
   char**
)
{
   if( !run(0, false, false) )
      return EXIT_FAILURE;

   if( !run(5, false, false) )
      return EXIT_FAILURE;

   if( !run(0, true, false) )
      return EXIT_FAILURE;

   if( !run(5, true, false) )
      return EXIT_FAILURE;

   if( !run(5, false, true) )
      return EXIT_FAILURE;

   std::cout << std::endl << "*** All tests passed" << std::endl;

   return EXIT_SUCCESS;
}
