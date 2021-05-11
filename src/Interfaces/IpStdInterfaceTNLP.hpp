// Copyright (C) 2004, 2010 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPSTDINTERFACETNLP_HPP__
#define __IPSTDINTERFACETNLP_HPP__

#include "IpUtils.hpp"
#include "IpTNLP.hpp"
#include "IpJournalist.hpp"
#include "IpException.hpp"
#include "IpStdCInterface.h"
#include "IpSmartPtr.hpp"

namespace Ipopt
{
/** Declare exception that is thrown when invalid NLP data is provided */
DECLARE_STD_EXCEPTION(INVALID_STDINTERFACE_NLP);

/** Implementation of a TNLP for the Standard C interface.
 *
 *  The standard C interface is exposed to the user as a single C
 *  function that is given problem dimension, starting points, and
 *  pointers for functions that evaluate objective function etc.
 */
class StdInterfaceTNLP : public TNLP
{
public:
   /**@name Constructors/Destructors */
   ///@{
   /** Constructor, given dimensions of problem, function pointers
    *  for evaluation callback functions, and starting points.
    *
    *  Note that the constructor does not make a copy of any of the Number
    *  arrays, i.e. it is up to the called to keep them around.
    */
   StdInterfaceTNLP(
      Index           n_var,
      const Number*   x_L,
      const Number*   x_U,
      Index           n_con,
      const Number*   g_L,
      const Number*   g_U,
      Index           nele_jac,
      Index           nele_hess,
      Index           index_style,
      const Number*   start_x,
      const Number*   start_lam,
      const Number*   start_z_L,
      const Number*   start_z_U,
      Eval_F_CB       eval_f,
      Eval_G_CB       eval_g,
      Eval_Grad_F_CB  eval_grad_f,
      Eval_Jac_G_CB   eval_jac_g,
      Eval_H_CB       eval_h,
      Intermediate_CB intermediate_cb,
      Number*         x_sol,
      Number*         z_L_sol,
      Number*         z_U_sol,
      Number*         g_sol,
      Number*         lam_sol,
      Number*         obj_sol,
      UserDataPtr     user_data,
      Number          obj_scaling = 1,
      const Number*   x_scaling = NULL,
      const Number*   g_scaling = NULL
   );

   /** Default destructor */
   virtual ~StdInterfaceTNLP();
   ///@}

   /**@name Methods to gather information about the NLP.
    *
    * These methods are overloaded from TNLP. See TNLP for their more detailed documentation.
    */
   ///@{
   virtual bool get_nlp_info(
      Index&          n,
      Index&          m,
      Index&          nnz_jac_g,
      Index&          nnz_h_lag,
      IndexStyleEnum& index_style
   );

   virtual bool get_bounds_info(
      Index   n,
      Number* x_l,
      Number* x_u,
      Index   m,
      Number* g_l,
      Number* g_u
   );

   virtual bool get_scaling_parameters(
      Number& obj_scaling,
      bool&   use_x_scaling,
      Index   n,
      Number* x_scaling,
      bool&   use_g_scaling,
      Index   m,
      Number* g_scaling
   );

   virtual bool get_starting_point(
      Index   n,
      bool    init_x,
      Number* x,
      bool    init_z,
      Number* z_L,
      Number* z_U,
      Index   m,
      bool    init_lambda,
      Number* lambda
   );

   virtual bool eval_f(
      Index         n,
      const Number* x,
      bool          new_x,
      Number&       obj_value
   );

   virtual bool eval_grad_f(
      Index         n,
      const Number* x,
      bool          new_x,
      Number*       grad_f
   );

   virtual bool eval_g(
      Index         n,
      const Number* x,
      bool          new_x,
      Index         m,
      Number*       g
   );

   virtual bool eval_jac_g(
      Index         n,
      const Number* x,
      bool          new_x,
      Index         m,
      Index         nele_jac,
      Index*        iRow,
      Index*        jCol,
      Number*       values
   );

   virtual bool eval_h(
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
      Number*       values
   );

   virtual bool intermediate_callback(
      AlgorithmMode              mode,
      Index                      iter,
      Number                     obj_value,
      Number                     inf_pr,
      Number                     inf_du,
      Number                     mu,
      Number                     d_norm,
      Number                     regularization_size,
      Number                     alpha_du,
      Number                     alpha_pr,
      Index                      ls_trials,
      const IpoptData*           ip_data,
      IpoptCalculatedQuantities* ip_cq
   );
   ///@}

   /** @name Solution Methods */
   ///@{
   virtual void finalize_solution(
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
      IpoptCalculatedQuantities* ip_cq
   );
   ///@}

   /// get_curr_iterate() to be called by GetIpoptCurrentIterate()
   /// @since 3.14.0
   bool get_curr_iterate(
      Bool                       scaled,
      Index                      n,
      Number*                    x,
      Number*                    z_L,
      Number*                    z_U,
      Index                      m,
      Number*                    g,
      Number*                    lambda
   ) const
   {
      return TNLP::get_curr_iterate(ip_data_, ip_cq_, scaled, n, x, z_L, z_U, m, g, lambda);
   }

   /// get_curr_violations() to be called by GetIpoptCurrentViolations()
   /// @since 3.14.0
   bool get_curr_violations(
      bool                       scaled,
      Index                      n,
      Number*                    x_L_violation,
      Number*                    x_U_violation,
      Number*                    compl_x_L,
      Number*                    compl_x_U,
      Number*                    grad_lag_x,
      Index                      m,
      Number*                    nlp_constraint_violation,
      Number*                    compl_g
   ) const
   {
      return TNLP::get_curr_violations(ip_data_, ip_cq_, scaled, n, x_L_violation, x_U_violation, compl_x_L, compl_x_U, grad_lag_x, m, nlp_constraint_violation, compl_g);
   }

private:
   /** Journalist */
   SmartPtr<const Journalist> jnlst_;

   /** @name Information about the problem */
   ///@{
   /** Number of variables */
   const Index n_var_;
   /** Number of constraints */
   const Index n_con_;
   /** Pointer to Number array containing lower bounds for variables */
   const Number* x_L_;
   /** Pointer to Number array containing upper bounds for variables */
   const Number* x_U_;
   /** Pointer to Number array containing lower bounds for constraints */
   const Number* g_L_;
   /** Pointer to Number array containing upper bounds for constraints */
   const Number* g_U_;
   /** Number of non-zero elements in the constraint Jacobian */
   const Index nele_jac_;
   /** Number of non-zero elements in the Hessian */
   const Index nele_hess_;
   /** Starting value of the iRow and jCol parameters for matrices */
   const Index index_style_;
   /** Pointer to Number array containing starting point for variables */
   const Number* start_x_;
   /** Pointer to Number array containing starting values for constraint multipliers */
   const Number* start_lam_;
   /** Pointer to Number array containing starting values for lower bound multipliers */
   const Number* start_z_L_;
   /** Pointer to Number array containing starting values for upper bound multipliers */
   const Number* start_z_U_;
   /** Pointer to callback function evaluating value of objective function */
   Eval_F_CB eval_f_;
   /** Pointer to callback function evaluating value of constraints */
   Eval_G_CB eval_g_;
   /** Pointer to callback function evaluating gradient of objective function */
   Eval_Grad_F_CB eval_grad_f_;
   /** Pointer to callback function evaluating Jacobian of constraints */
   Eval_Jac_G_CB eval_jac_g_;
   /** Pointer to callback function evaluating Hessian of Lagrangian */
   Eval_H_CB eval_h_;
   /** Pointer to intermediate callback function giving control to user */
   Intermediate_CB intermediate_cb_;
   /** Pointer to user data */
   UserDataPtr user_data_;
   /** Objective scaling factor */
   Number obj_scaling_;
   /** Scaling factors for variables (if not NULL) */
   const Number* x_scaling_;
   /** Scaling factors for constraints (if not NULL) */
   const Number* g_scaling_;
   ///@}

   /** A non-const copy of x - this is kept up-to-date in apply_new_x */
   Number* non_const_x_;

   /** @name Pointers to the user provided vectors for solution */
   ///@{
   Number* x_sol_;
   Number* z_L_sol_;
   Number* z_U_sol_;
   Number* g_sol_;
   Number* lambda_sol_;
   Number* obj_sol_;
   ///@}

   /** @name Temporary pointers to IpoptData and IpoptCalculatedQuantities
    *
    * For implementation of GetIpoptCurrentIterate() and GetIpoptCurrentViolations() (without API change).
    */
   const IpoptData*           ip_data_;
   IpoptCalculatedQuantities* ip_cq_;

   /** Update the internal state if the x value changes */
   void apply_new_x(
      bool          new_x,
      Index         n,
      const Number* x
   );

   /**@name Default Compiler Generated Methods
    * (Hidden to avoid implicit creation/calling).
    * These methods are not implemented and
    * we do not want the compiler to implement
    * them for us, so we declare them private
    * and do not define them. This ensures that
    * they will not be implicitly created/called. */
   ///@{
   /** Default Constructor */
   StdInterfaceTNLP();

   /** Copy Constructor */
   StdInterfaceTNLP(
      const StdInterfaceTNLP&
   );

   /** Default Assignment Operator */
   void operator=(
      const StdInterfaceTNLP&
   );
   ///@}

};

} // namespace Ipopt

#endif
