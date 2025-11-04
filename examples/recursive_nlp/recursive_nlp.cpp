// Copyright (C) 2020 COIN-OR Foundation
// All Rights Reserved.
// This code is published under the Eclipse Public License 2.0.
//
// Author:  Brad Bell

/** @file recursive_nlp.cpp
 *
 * This solve a NLP where the objective function is defined as a solution point of a NLP.
 * Both the outer and the inner problem are solved by Ipopt.
 * (For simplicity, the inner problem is very similar to the outer problem.)
 *
 * Inner problem:  minimize exp [ 0.5 * (x - a)^2    + 0.5 * x^2 ] w.r.t x
 * Define y(a):    arg_min  exp [ 0.5 * (x - a)^2    + 0.5 * x^2 ] w.r.t x
 * Outer problem:  minimize exp [ 0.5 * (y(a) - a)^2 + 0.5 * y(a)^2 ] w.r.t a
 */

#include "IpIpoptApplication.hpp"
#include "IpTNLP.hpp"

#include <iostream>
#include <cassert>
#include <cmath>

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

using namespace Ipopt;

class recursive_nlp: public TNLP
{
private:
   const bool   inner_;
   const Number a_;
   Number       arg_min_;
   Number       y_;

public:
   // arg_min
   Number arg_min(void) const
   {
      return arg_min_;
   }

   // constructor for the inner problem
   recursive_nlp(const Number& a)
      : inner_(true),
        a_(a),
        arg_min_(0.0),
        y_(0.0)
   { }

   // constructor for the outer problem
   recursive_nlp(void)
      : inner_(false),
        a_(0.0),
        arg_min_(0.0),
        y_(0.0)
   { }

   // default destructor
   virtual ~recursive_nlp() { }

   // get_nlp_info
   bool get_nlp_info(
      Index&          n,
      Index&          m,
      Index&          nnz_jac_g,
      Index&          nnz_h_lag,
      IndexStyleEnum& index_style)
   {
      n = 1;
      m = 0;
      nnz_jac_g = 0;
      nnz_h_lag = 1;
      index_style = C_STYLE;

      return true;
   }

   // get_bounds_info
   bool get_bounds_info(
      Index   n,
      Number* x_l,
      Number* x_u,
      Index   m,
      Number* g_l,
      Number* g_u)
   {
      assert(n == 1 && m == 0);

      /* default Ipopt value for "infinity" is 1e20 */
      x_l[0] = -1e20;
      x_u[0] = +1e20;

      return true;
   }

   // get_starting_point
   bool get_starting_point(
      Index   n,
      bool    init_x,
      Number* x,
      bool    init_z,
      Number* z_L,
      Number* z_U,
      Index   m,
      bool    init_lambda,
      Number* lambda)
   {
      assert(n == 1 && m == 0 && init_x);

      // x[0] == 0 is solution for outer problem
      x[0] = 2.0;

      return true;
   }

   // eval_f
   bool eval_f(
      Index         n,
      const Number* x,
      bool          new_x,
      Number&       obj_value)
   {
      assert(n == 1);

      if( inner_ )
      {
         Number arg = 0.5 * (x[0] - a_) * (x[0] - a_) + 0.5 * x[0] * x[0];
         // avoid a floating-point overflow when arg is too large, return false instead
         if( arg >= std::log(std::numeric_limits<Number>::max()) )
         {
            return false;
         }
         obj_value = std::exp(arg);
         return true;
      }

      if( new_x )
      {
         // solve the inner problem with corresponding to x[0]
         SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
         app->Options()->SetIntegerValue("print_level", J_STRONGWARNING);
         app->Options()->SetStringValue("hessian_approximation", "limited-memory");
         ApplicationReturnStatus status = app->Initialize();
         if( status != Solve_Succeeded )
         {
            return false;
         }

         SmartPtr<recursive_nlp> nlp = new recursive_nlp(x[0]);
         status = app->OptimizeTNLP(nlp);
         if( status != Solve_Succeeded )
         {
            return false;
         }

         // set y_ equal to the arg_min for the inner problem
         y_ = nlp->arg_min();
      }

      // evaluate object for the outer problem
      Number arg = 0.5 * (y_ - x[0]) * (y_ - x[0]) + 0.5 * y_ * y_;
      // avoid a floating-point overflow when arg is too large, return false instead
      if( arg >= std::log(std::numeric_limits<Number>::max()) )
      {
         return false;
      }
      obj_value = std::exp(arg);

      return true;
   }

   // eval_grad_f
   bool eval_grad_f(
      Index         n,
      const Number* x,
      bool          new_x,
      Number*       grad_f)
   {
      assert(n == 1);

      // evaluate objective, if new_x and outer problem also update y_
      Number obj_value;
      eval_f(n, x, new_x, obj_value);

      if( inner_ )
      {
         grad_f[0] = obj_value * ((x[0] - a_) + x[0]);
         return true;
      }

      // use fact that gradient of outer objective w.r.t y_ is zero
      grad_f[0] = obj_value * (x[0] - y_);

      return true;
   }

   // eval_g
   bool eval_g(
      Index         n,
      const Number* x,
      bool          new_x,
      Index         m,
      Number*       g)
   {
      assert(n == 1 && m == 0);

      return true;
   }

   // eval_jac_g
   bool eval_jac_g(
      Index         n,
      const Number* x,
      bool          new_x,
      Index         m,
      Index         nele_jac,
      Index*        iRow,
      Index*        jCol,
      Number*       values)
   {
      assert(n == 1 && m == 0 && nele_jac == 0);

      return true;
   }

   // eval_h
   bool eval_h(
      Index         n,
      const Number* x,
      bool          new_x,
      Number        obj_factor,
      Index         m,
      const Number* lambda,
      bool          new_lambda,
      Index         nele_hess,
      Index*        iRow,
      Index*        jCol,
      Number*       values)
   {
      // using an approximate Hessian
      return false;
   }

   // finalize_solution
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
      const IpoptData*           ip_data,
      IpoptCalculatedQuantities* ip_cq)
   {
      assert(n == 1 && m == 0);

      arg_min_ = x[0];
   }
};

int main(int argc, char** argv)
{
   SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
   ApplicationReturnStatus status = app->Initialize();
   assert(status == Solve_Succeeded);

   app->Options()->SetStringValue("derivative_test", "first-order");
   app->Options()->SetStringValue("hessian_approximation", "limited-memory");

   SmartPtr<recursive_nlp> nlp = new recursive_nlp();
   status = app->OptimizeTNLP(nlp);

   if( status == Solve_Succeeded )
   {
      std::cout << "nlp->arg_min() = " << nlp->arg_min() << "\n";
   }

   return EXIT_SUCCESS;
}
