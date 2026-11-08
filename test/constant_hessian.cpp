// This code is published under the Eclipse Public License.

#include "IpIpoptApplication.hpp"
#include "IpTNLP.hpp"

#include <cstdlib>
#include <iostream>

using namespace Ipopt;

// min coefficient * x^2 / 2 - x, subject to -1 <= x <= 1.
class ConstantHessianNLP: public TNLP
{
public:
   Number coefficient;
   Index hessian_evaluations;

   ConstantHessianNLP()
      : coefficient(2.0),
        hessian_evaluations(0)
   { }

   bool get_nlp_info(
      Index&          n,
      Index&          m,
      Index&          nnz_jac_g,
      Index&          nnz_h_lag,
      IndexStyleEnum& index_style
   )
   {
      n = 1;
      m = 0;
      nnz_jac_g = 0;
      nnz_h_lag = 1;
      index_style = C_STYLE;
      return true;
   }

   bool get_bounds_info(
      Index,
      Number* x_l,
      Number* x_u,
      Index,
      Number*,
      Number*
   )
   {
      x_l[0] = -1.0;
      x_u[0] = 1.0;
      return true;
   }

   bool get_starting_point(
      Index,
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
      {
         x[0] = 0.8;
      }
      return !init_z && !init_lambda;
   }

   bool eval_f(
      Index,
      const Number* x,
      bool,
      Number&       obj_value
   )
   {
      obj_value = coefficient * x[0] * x[0] / 2.0 - x[0];
      return true;
   }

   bool eval_grad_f(
      Index,
      const Number* x,
      bool,
      Number*       grad_f
   )
   {
      grad_f[0] = coefficient * x[0] - 1.0;
      return true;
   }

   bool eval_g(Index, const Number*, bool, Index, Number*)
   {
      return true;
   }

   bool eval_jac_g(Index, const Number*, bool, Index, Index, Index*, Index*, Number*)
   {
      return true;
   }

   bool eval_h(
      Index,
      const Number*,
      bool,
      Number  obj_factor,
      Index,
      const Number*,
      bool,
      Index,
      Index*  iRow,
      Index*  jCol,
      Number* values
   )
   {
      if( values == NULL )
      {
         iRow[0] = jCol[0] = 0;
      }
      else
      {
         ++hessian_evaluations;
         values[0] = obj_factor * coefficient;
      }
      return true;
   }

   void finalize_solution(
      SolverReturn,
      Index,
      const Number*,
      const Number*,
      const Number*,
      Index,
      const Number*,
      const Number*,
      Number,
      const IpoptData*,
      IpoptCalculatedQuantities*
   )
   { }
};

int main()
{
   SmartPtr<ConstantHessianNLP> nlp = new ConstantHessianNLP();
   SmartPtr<IpoptApplication> app = new IpoptApplication();
   app->Options()->SetStringValue("hessian_approximation", "exact");
   app->Options()->SetStringValue("hessian_constant", "yes");
   app->Options()->SetStringValue("nlp_scaling_method", "none");
   app->Options()->SetIntegerValue("print_level", 0);
   app->Options()->SetNumericValue("tol", 1e-5);

   if( app->Initialize() != Solve_Succeeded )
   {
      return EXIT_FAILURE;
   }

   for( Index solve = 0; solve < 2; ++solve )
   {
      ApplicationReturnStatus status;
      if( solve == 0 )
      {
         status = app->OptimizeTNLP(nlp);
      }
      else
      {
         // Change the numerical Hessian, keeping the same TNLP and sparsity.
         nlp->coefficient = 10.0;
         nlp->hessian_evaluations = 0;
         app->Options()->SetStringValue("warm_start_same_structure", "yes");
         status = app->ReOptimizeTNLP(nlp);
      }

      // Refresh between solves, but retain constant-Hessian reuse within each solve.
      if( nlp->hessian_evaluations != 1 )
      {
         std::cerr << "Solve " << solve + 1 << ": expected one Hessian evaluation, got "
                   << nlp->hessian_evaluations << std::endl;
         return EXIT_FAILURE;
      }
      if( status != Solve_Succeeded )
      {
         std::cerr << "Solve " << solve + 1 << " returned " << status << std::endl;
         return EXIT_FAILURE;
      }
   }

   return EXIT_SUCCESS;
}
