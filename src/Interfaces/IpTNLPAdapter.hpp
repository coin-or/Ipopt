// Copyright (C) 2004, 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPTNLPADAPTER_HPP__
#define __IPTNLPADAPTER_HPP__

#include "IpNLP.hpp"
#include "IpTNLP.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpExpansionMatrix.hpp"
#include <list>

namespace Ipopt
{

// forward declarations
class ExpansionMatrixSpace;
class IteratesVector;
class TDependencyDetector;

/** This class adapts the TNLP interface so it looks like an NLP interface.
 *
 *  This is an Adapter class (Design Patterns) that converts a TNLP to an NLP.
 *  This allows users to write to the "more convenient" TNLP interface.
 *
 *  Given a TNLP
 *  \f{eqnarray*}
 *     \mathrm{min}  && f(x), \\
 *     \mathrm{s.t.} && g_L \leq g(x) \leq g_U, &\qquad \lambda\\
 *                   && x_L \leq  x   \leq x_U, &\qquad z_L, z_U
 *  \f}
 *  let \f$E = \{i : g_{L,i} = g_{U,i}\}\f$ and \f$I = \{i : g_{L,i} \neq g_{U,i}\}\f$
 *  be the indices of equality and inequality constraints, respectively.
 *  The dual variables for the constraints are \f$\lambda\f$.
 *  The dual variables for the variable bounds are \f$z_L\f$ and \f$z_U\f$.
 *
 *  A TNLPAdapter represents the problem
 *  \f{eqnarray*}
 *     \mathrm{min}  && f(x), \\
 *     \mathrm{s.t.} && c(x) = 0,               &\qquad y_c\\
 *                   && d_L \leq d(x) \leq d_U, &\qquad y_d \\
 *                   && x_L \leq  x \leq x_U,   &\qquad z_L, z_U
 *  \f}
 *  where \f$c(x) = g_E(x) - g_{L,E}\f$, i.e., corresponding to equality constraints of TNLP, and
 *  \f$d(x) = g_I(x)\f$, \f$d_L = g_{L,I}\f$, \f$d_U = g_{U,I}\f$, i.e., corresponding to inequality constraints of TNLP.
 *
 *  The dual variables for the constraints are \f$y_c\f$ and \f$y_d\f$.
 *  The dual variables for the bounds on slack and original variables are \f$s_L\f$, \f$s_U\f$, \f$z_L\f$, \f$z_U\f$.
 *
 *  Internally, %Ipopt reformulates \f$d_L \leq d(x) \leq d_U\f$ as
 *  \f{eqnarray*}
 *                   && d(x) - s = 0        &\qquad y_d\\
 *                   && d_L \leq s \leq d_U &\qquad v_L, v_U \\
 *  \f}
 *
 *  If fixed variables are present (\f$x_{L,i} = x_{U,i}\f$) in the TNLP and
 *  fixed_variable_treatment is set to make_parameter, then these variables are not made visible
 *  to Ipopt internally.
 *  If fixed_variable_treatment is set to make_constraint, then their bounds relaxed and equality
 *  constraints \f$x_i - x_{L,i} = 0\f$ are added to the end of \f$c(x) = 0\f$.
 */
class IPOPTLIB_EXPORT TNLPAdapter : public NLP
{
public:
   /**@name Constructors/Destructors */
   ///@{
   /** Default constructor */
   TNLPAdapter(
      const SmartPtr<TNLP>             tnlp,
      const SmartPtr<const Journalist> jnlst = NULL
   );

   /** Default destructor */
   virtual ~TNLPAdapter();
   ///@}

   /**@name Exceptions */
   ///@{
   DECLARE_STD_EXCEPTION(INVALID_TNLP);
   DECLARE_STD_EXCEPTION(ERROR_IN_TNLP_DERIVATIVE_TEST);
   ///@}

   /** @name TNLPAdapter Initialization. */
   ///@{
   virtual bool ProcessOptions(
      const OptionsList& options,
      const std::string& prefix
   );

   /** Method for creating the derived vector / matrix types.
    *
    *  Do not delete these.
    */
   virtual bool GetSpaces(
      SmartPtr<const VectorSpace>&    x_space,
      SmartPtr<const VectorSpace>&    c_space,
      SmartPtr<const VectorSpace>&    d_space,
      SmartPtr<const VectorSpace>&    x_l_space,
      SmartPtr<const MatrixSpace>&    px_l_space,
      SmartPtr<const VectorSpace>&    x_u_space,
      SmartPtr<const MatrixSpace>&    px_u_space,
      SmartPtr<const VectorSpace>&    d_l_space,
      SmartPtr<const MatrixSpace>&    pd_l_space,
      SmartPtr<const VectorSpace>&    d_u_space,
      SmartPtr<const MatrixSpace>&    pd_u_space,
      SmartPtr<const MatrixSpace>&    Jac_c_space,
      SmartPtr<const MatrixSpace>&    Jac_d_space,
      SmartPtr<const SymMatrixSpace>& Hess_lagrangian_space
   );

   /** Method for obtaining the bounds information. */
   virtual bool GetBoundsInformation(
      const Matrix& Px_L,
      Vector&       x_L,
      const Matrix& Px_U,
      Vector&       x_U,
      const Matrix& Pd_L,
      Vector&       d_L,
      const Matrix& Pd_U,
      Vector&       d_U
   );

   /** Method for obtaining the starting point for all the iterates. */
   virtual bool GetStartingPoint(
      SmartPtr<Vector> x,
      bool             need_x,
      SmartPtr<Vector> y_c,
      bool             need_y_c,
      SmartPtr<Vector> y_d,
      bool             need_y_d,
      SmartPtr<Vector> z_L,
      bool             need_z_L,
      SmartPtr<Vector> z_U,
      bool             need_z_U
   );

   /** Method for obtaining an entire iterate as a warmstart point.
    *
    * The incoming IteratesVector has to be filled.
    */
   virtual bool GetWarmStartIterate(
      IteratesVector& warm_start_iterate
   );
   ///@}

   /** @name TNLPAdapter evaluation routines. */
   ///@{
   virtual bool Eval_f(
      const Vector& x,
      Number&       f
   );

   virtual bool Eval_grad_f(
      const Vector& x,
      Vector&       g_f
   );

   virtual bool Eval_c(
      const Vector& x,
      Vector&       c
   );

   virtual bool Eval_jac_c(
      const Vector& x,
      Matrix&       jac_c
   );

   virtual bool Eval_d(
      const Vector& x,
      Vector&       d
   );

   virtual bool Eval_jac_d(
      const Vector& x,
      Matrix&       jac_d
   );

   virtual bool Eval_h(
      const Vector& x,
      Number        obj_factor,
      const Vector& yc,
      const Vector& yd,
      SymMatrix&    h
   );

   virtual void GetScalingParameters(
      const SmartPtr<const VectorSpace> x_space,
      const SmartPtr<const VectorSpace> c_space,
      const SmartPtr<const VectorSpace> d_space,
      Number&                           obj_scaling,
      SmartPtr<Vector>&                 x_scaling,
      SmartPtr<Vector>&                 c_scaling,
      SmartPtr<Vector>&                 d_scaling
   ) const;
   ///@}

   /** @name Solution Reporting Methods */
   ///@{
   virtual void FinalizeSolution(
      SolverReturn               status,
      const Vector&              x,
      const Vector&              z_L,
      const Vector&              z_U,
      const Vector&              c,
      const Vector&              d,
      const Vector&              y_c,
      const Vector&              y_d,
      Number                     obj_value,
      const IpoptData*           ip_data,
      IpoptCalculatedQuantities* ip_cq
   );

   virtual bool IntermediateCallBack(
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

   /** Method returning information on quasi-Newton approximation. */
   virtual void GetQuasiNewtonApproximationSpaces(
      SmartPtr<VectorSpace>& approx_space,
      SmartPtr<Matrix>&      P_approx
   );

   /** Enum for treatment of fixed variables option */
   enum FixedVariableTreatmentEnum
   {
      MAKE_PARAMETER = 0,
      MAKE_PARAMETER_NODUAL,  ///< @since 3.14.0
      MAKE_CONSTRAINT,
      RELAX_BOUNDS
   };

   /** Enum for specifying which derivative test is to be performed. */
   enum DerivativeTestEnum
   {
      NO_TEST = 0,
      FIRST_ORDER_TEST,
      SECOND_ORDER_TEST,
      ONLY_SECOND_ORDER_TEST
   };

   /** Enum for specifying technique for computing Jacobian */
   enum JacobianApproxEnum
   {
      JAC_EXACT = 0,
      JAC_FINDIFF_VALUES
   };

   /** Enum for specifying technique for computing objective Gradient */
   enum GradientApproxEnum
   {
      OBJGRAD_EXACT = 0,
      OBJGRAD_FINDIFF_VALUES
   };

   /** Method for performing the derivative test */
   bool CheckDerivatives(
      DerivativeTestEnum deriv_test,
      Index              deriv_test_start_index
   );

   static void RegisterOptions(
      SmartPtr<RegisteredOptions> roptions
   );

   /** Accessor method for the underlying TNLP. */
   SmartPtr<TNLP> tnlp() const
   {
      return tnlp_;
   }

   /** @name Methods for translating data for IpoptNLP into the TNLP data.
    *
    *  These methods are used to obtain the current (or final)
    *  data for the TNLP formulation from the IpoptNLP structure.
    */
   ///@{
   /** Sort the primal variables, and add the fixed values in x_orig */
   void ResortX(
      const Vector& x,                   /**< internal values for primal variables x */
      Number*       x_orig,              /**< vector to fill with values from x */
      bool          usefixedvals = true  /**< whether to use stored variable fixings for fixed variables (true), or zero (false) @since 3.14.0 */
   );

   /** Sort constraint values
    *
    * To Ipopt, the equality constraints are presented with right hand side zero.
    * Specifying correctrhs=true corrects for the original right hand side.
    */
   void ResortG(
      const Vector& c,                   /**< internal activity for equality constraints */
      const Vector& d,                   /**< internal activity for inequality constraints */
      Number*       g_orig,              /**< vector to fill with values from c and d */
      bool          correctrhs = false   /**< whether to add unscaled rhs-values for constraints that internally correspond to c(x)=0 @since 3.14.0 */
   );

   /** Provides values for lower and upper bounds on variables for given Ipopt-internal vectors.
    *
    * Similar to ResortX, but does so for two arrays and does not set any values for fixed variables.
    * @since 3.14.0
    */
   void ResortBounds(
      const Vector& x_L,              /**< internal values for lower bounds on x */
      Number*       x_L_orig,         /**< vector to fill with values from x_L */
      const Vector& x_U,              /**< internal values for upper bounds on x */
      Number*       x_U_orig          /**< vector to fill with values from x_U */
   );

   /** Provides values for lower and upper bounds on variables for given Ipopt-internal vectors.
    *
    * Similar to ResortX, but does so for two arrays and does not set any values for fixed variables.
    *
    * @deprecated Use ResortBounds() instead.
    */
   IPOPT_DEPRECATED
   void ResortBnds(
      const Vector& x_L,              /**< internal values for lower bounds on x */
      Number*       x_L_orig,         /**< vector to fill with values from x_L */
      const Vector& x_U,              /**< internal values for upper bounds on x */
      Number*       x_U_orig,         /**< vector to fill with values from x_U */
      bool          clearorig = true  /**< whether to initialize complete x_L_orig and x_U_orig to 0.0 before setting values for non-fixed variables - ignored */
   )
   {
      ResortBounds(x_L, x_L_orig, x_U, x_U_orig);
      (void) clearorig;
   }

   /** Provides values for dual multipliers on lower and upper bounds on variables for given Ipopt-internal vectors.
    *
    * Similar to ResortBounds, but also provides dual values for fixed variables if fixed_variable_treatment is set to make_constraint or make_parameter.
    *
    * @attention If there are fixed variables and fixed_variable_treatment is make_parameter (the default),
    *   then the Gradient of f(x) and the Jacobian of g(x) may be reevaluated here (that's why the function needs x).
    *   Further, in this setting, this function only provides correct bound multipliers for fixed variables
    *   if x, y_c, and y_d correspond to the unscaled problem.
    *
    * @return True, if bound multipliers could be assigned. False if there was an evaluation error when calculating bound multipliers for fixed variables.
    * @since 3.14.0
    */
   bool ResortBoundMultipliers(
      const Vector& x,                /**< internal values for primal variables x */
      const Vector& y_c,              /**< internal values for equality constraint multipliers */
      const Vector& y_d,              /**< internal values for inequality constraint multipliers */
      const Vector& z_L,              /**< internal values for lower bound multipliers */
      Number*       z_L_orig,         /**< vector to fill with values from z_L */
      const Vector& z_U,              /**< internal values for upper bound multipliers */
      Number*       z_U_orig          /**< vector to fill with values from z_U */
   );

   /** Get number of variables and number of constraints in TNLP
    * @since 3.14.0
    */
   void GetFullDimensions(
      Index& n,  /**< storage for full dimension of x (fixed + non-fixed) */
      Index& m   /**< storage for full dimension of g (c + d) */
   ) const
   {
      n = n_full_x_;
      m = n_full_g_;
   }

   /** Get number and indices of fixed variables
    * @since 3.14.0
    */
   void GetFixedVariables(
      Index&  n_x_fixed,     /**< storage for number of fixed variables in TNLP */
      Index*& x_fixed_map,   /**< storage for pointer to array that holds indices of fixed variables (has length n_fixed_x, can be NULL if n_fixed_x=0) */
      FixedVariableTreatmentEnum& fixed_variable_treatment  /**< treatment for fixed variables as used by TNLP */
   ) const
   {
      n_x_fixed = n_x_fixed_;
      x_fixed_map = x_fixed_map_;
      fixed_variable_treatment = fixed_variable_treatment_;
   }

   /** Get mappings between TNLP indices and Ipopt internal indices for variables and constraints.
    *
    * See the various ResortXyz functions on usage.
    *
    * \attention P_x_full_x is set to NULL if there are no fixed variables or fixed_variable_treatment is not make_parameter
    * @since 3.14.0
    */
   void GetPermutationMatrices(
      SmartPtr<const ExpansionMatrix>& P_x_full_x, /**< map from TNLP variable indices to Ipopt internal indices, filtered variables get mapped to -1 */
      SmartPtr<const ExpansionMatrix>& P_x_x_L,    /**< map from indices on lower bounds on x to Ipopt internal indices for x */
      SmartPtr<const ExpansionMatrix>& P_x_x_U,    /**< map from indices on upper bounds on x to Ipopt internal indices for x */
      SmartPtr<const ExpansionMatrix>& P_c_g,      /**< map from indices on equality constraints (c(x)=0) into TNLP constraint indices (g_l <= g(x) <= g_u) */
      SmartPtr<const ExpansionMatrix>& P_d_g       /**< map from indices on inequality constraints (d(x)-s=0) into TNLP constraint indices (g_l <= g(x) <= g_u) */
   ) const
   {
      P_x_full_x = ConstPtr(P_x_full_x_);
      P_x_x_L = ConstPtr(P_x_x_L_);
      P_x_x_U = ConstPtr(P_x_x_U_);
      P_c_g = ConstPtr(P_c_g_);
      P_d_g = ConstPtr(P_d_g_);
   }

   /** Get right-hand-sides that are added into c(x)
    * @since 3.14.0
    */
   const Number* GetC_Rhs() const
   {
      return c_rhs_;
   }

   ///@}

private:
   /**@name Default Compiler Generated Methods
    * (Hidden to avoid implicit creation/calling).
    *
    * These methods are not implemented and
    * we do not want the compiler to implement
    * them for us, so we declare them private
    * and do not define them. This ensures that
    * they will not be implicitly created/called.
    */
   ///@{
   /** Copy Constructor */
   TNLPAdapter(
      const TNLPAdapter&
   );

   /** Default Assignment Operator */
   void operator=(
      const TNLPAdapter&
   );
   ///@}

   /** @name Method implementing the detection of linearly dependent equality constraints */
   bool DetermineDependentConstraints(
      Index             n_x_var,
      const Index*      x_not_fixed_map,
      const Number*     x_l,
      const Number*     x_u,
      const Number*     g_l,
      const Number*     g_u,
      Index             n_c,
      const Index*      c_map,
      std::list<Index>& c_deps);

   /** Pointer to the TNLP class (class specific to Number* vectors and triplet matrices) */
   SmartPtr<TNLP> tnlp_;

   /** Journalist */
   SmartPtr<const Journalist> jnlst_;

   /** Object that can be used to detect linearly dependent rows in the equality constraint Jacobian */
   SmartPtr<TDependencyDetector> dependency_detector_;

   /**@name Algorithmic parameters */
   ///@{
   /** Value for a lower bound that denotes -infinity */
   Number nlp_lower_bound_inf_;
   /** Value for a upper bound that denotes infinity */
   Number nlp_upper_bound_inf_;
   /** Flag indicating how fixed variables should be handled */
   FixedVariableTreatmentEnum fixed_variable_treatment_;
   /** Determines relaxation of fixing bound for RELAX_BOUNDS. */
   Number bound_relax_factor_;
   /** Maximal slack for one-sidedly bounded variables.  If a
    *  variable has only one bound, say a lower bound xL, then an
    *  upper bound xL + max_onesided_bound_slack_.  If this value is
    *  zero, no upper bound is added. */
   /* Took this out:  Number max_onesided_bound_slack_; */
   /** Whether and which derivative test should be performed at starting point */
   DerivativeTestEnum derivative_test_;
   /** Size of the perturbation for the derivative test */
   Number derivative_test_perturbation_;
   /** Relative threshold for marking deviation from finite difference test */
   Number derivative_test_tol_;
   /** Flag indicating if all test values should be printed, or only those violating the threshold. */
   bool derivative_test_print_all_;
   /** Index of first quantity to be checked. */
   Index derivative_test_first_index_;
   /** Flag indicating whether the TNLP with identical structure has already been solved before. */
   bool warm_start_same_structure_;
   /** Flag indicating what Hessian information is to be used. */
   HessianApproximationType hessian_approximation_;
   /** Number of linear variables. */
   Index num_linear_variables_;
   /** Flag indicating how Jacobian is computed. */
   JacobianApproxEnum jacobian_approximation_;
   /** Flag indicating how objective Gradient is computed. */
   GradientApproxEnum gradient_approximation_;
   /** Size of the perturbation for the derivative approximation */
   Number findiff_perturbation_;
   /** Maximal perturbation of the initial point */
   Number point_perturbation_radius_;
   /** Flag indicating if rhs should be considered during dependency detection */
   bool dependency_detection_with_rhs_;

   /** Overall convergence tolerance */
   Number tol_;
   ///@}

   /**@name Problem Size Data */
   ///@{
   /** full dimension of x (fixed + non-fixed) */
   Index n_full_x_;
   /** full dimension of g (c + d) */
   Index n_full_g_;
   /** non-zeros of the jacobian of c */
   Index nz_jac_c_;
   /** non-zeros of the jacobian of c without added constraints for fixed variables. */
   Index nz_jac_c_no_extra_;
   /** non-zeros of the jacobian of d */
   Index nz_jac_d_;
   /** number of non-zeros in full-size Jacobian of g */
   Index nz_full_jac_g_;
   /** number of non-zeros in full-size Hessian */
   Index nz_full_h_;
   /** number of non-zeros in the non-fixed-size Hessian */
   Index nz_h_;
   /** Number of fixed variables */
   Index n_x_fixed_;
   ///@}

   /** Numbering style of variables and constraints */
   TNLP::IndexStyleEnum index_style_;

   /** @name Local copy of spaces (for warm start) */
   ///@{
   SmartPtr<const VectorSpace> x_space_;
   SmartPtr<const VectorSpace> c_space_;
   SmartPtr<const VectorSpace> d_space_;
   SmartPtr<const VectorSpace> x_l_space_;
   SmartPtr<const MatrixSpace> px_l_space_;
   SmartPtr<const VectorSpace> x_u_space_;
   SmartPtr<const MatrixSpace> px_u_space_;
   SmartPtr<const VectorSpace> d_l_space_;
   SmartPtr<const MatrixSpace> pd_l_space_;
   SmartPtr<const VectorSpace> d_u_space_;
   SmartPtr<const MatrixSpace> pd_u_space_;
   SmartPtr<const MatrixSpace> Jac_c_space_;
   SmartPtr<const MatrixSpace> Jac_d_space_;
   SmartPtr<const SymMatrixSpace> Hess_lagrangian_space_;
   ///@}

   /**@name Local Copy of the Data */
   ///@{
   Number* full_x_; /** copy of the full x vector (fixed & non-fixed) */
   Number* full_lambda_; /** copy of lambda (yc & yd) */
   Number* full_g_; /** copy of g (c & d) */
   Number* jac_g_; /** the values for the full jacobian of g */
   Number* c_rhs_; /** the rhs values of c */
   ///@}

   /**@name Tags for deciding when to update internal copies of vectors */
   ///@{
   TaggedObject::Tag x_tag_for_iterates_;
   TaggedObject::Tag y_c_tag_for_iterates_;
   TaggedObject::Tag y_d_tag_for_iterates_;
   TaggedObject::Tag x_tag_for_g_;
   TaggedObject::Tag x_tag_for_jac_g_;
   ///@}

   /**@name Methods to update the values in the local copies of vectors */
   ///@{
   bool update_local_x(const Vector& x);
   bool update_local_lambda(const Vector& y_c, const Vector& y_d);
   ///@}

   /**@name Internal routines for evaluating g and jac_g.
    *
    * Values stored since they are used in both c and d routines.
    */
   ///@{
   bool internal_eval_g(bool new_x);
   bool internal_eval_jac_g(bool new_x);
   ///@}

   /** @name Internal methods for dealing with finite difference approximation */
   ///@{
   /** Initialize sparsity structure for finite difference Jacobian */
   void initialize_findiff_jac(const Index* iRow, const Index* jCol);
   ///@}

   /**@name Internal Permutation Spaces and matrices
    */
   ///@{
   /** Expansion from fixed x (ipopt) to full x */
   SmartPtr<ExpansionMatrix> P_x_full_x_;
   SmartPtr<ExpansionMatrixSpace> P_x_full_x_space_;

   /** Expansion from fixed x_L (ipopt) to full x */
   SmartPtr<ExpansionMatrix> P_x_x_L_;
   SmartPtr<ExpansionMatrixSpace> P_x_x_L_space_;

   /** Expansion from fixed x_U (ipopt) to full x */
   SmartPtr<ExpansionMatrix> P_x_x_U_;
   SmartPtr<ExpansionMatrixSpace> P_x_x_U_space_;

   /** Expansion from c only (ipopt) to full ampl c */
   SmartPtr<ExpansionMatrixSpace> P_c_g_space_;
   SmartPtr<ExpansionMatrix> P_c_g_;

   /** Expansion from d only (ipopt) to full ampl d */
   SmartPtr<ExpansionMatrixSpace> P_d_g_space_;
   SmartPtr<ExpansionMatrix> P_d_g_;

   Index* jac_idx_map_;
   Index* h_idx_map_;

   /** Position of fixed variables. This is required for a warm start */
   Index* x_fixed_map_;

   /** Index mapping of Jacobian w.r.t. fixed variables. */
   std::vector<Index> jac_fixed_idx_map_;
   std::vector<Index> jac_fixed_iRow_;
   std::vector<Index> jac_fixed_jCol_;
   ///@}

   /** @name Data for finite difference approximations of derivatives */
   ///@{
   /** Number of unique nonzeros in constraint Jacobian */
   Index findiff_jac_nnz_;
   /** Start position for nonzero indices in ja for each column of Jacobian */
   Index* findiff_jac_ia_;
   /** Ordered by columns, for each column the row indices in Jacobian */
   Index* findiff_jac_ja_;
   /** Position of entry in original triplet matrix */
   Index* findiff_jac_postriplet_;
   /** Copy of the lower bounds */
   Number* findiff_x_l_;
   /** Copy of the upper bounds */
   Number* findiff_x_u_;
   ///@}
};

} // namespace Ipopt

#endif
