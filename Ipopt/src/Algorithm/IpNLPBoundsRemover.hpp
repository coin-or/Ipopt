// Copyright (C) 2008, 2010 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: IpNLP.hpp 949 2007-03-27 00:41:26Z andreasw $
//
// Authors:  Andreas Waechter           IBM    2008-08-25

#ifndef __IPNLPBOUNDSREMOVER_HPP__
#define __IPNLPBOUNDSREMOVER_HPP__

#include "IpNLP.hpp"

namespace Ipopt
{
  /** This is an adaper for an NLP that converts variable bound
   *  constraints to inequality constraints.  This is necessary for
   *  the version of Ipopt that uses iterative linear solvers.  At
   *  this point, none of the original inequality constraints is
   *  allowed to have both lower and upper bounds.  The NLP visible to
   *  Ipopt via this adapter will not have any bounds on variables,
   *  but have equivalent inequality constraints.
   */
  class NLPBoundsRemover : public NLP
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** The constructor is given the NLP of which the bounds are to be
     *  replaced by inequality constriants.  */
    NLPBoundsRemover(NLP& nlp, bool allow_twosided_inequalities = false);

    /** Default destructor */
    virtual ~NLPBoundsRemover()
    {}
    //@}

    /** @name NLP Initialization (overload in
     *  derived classes).*/
    //@{
    /** Overload if you want the chance to process options or parameters that
     *  may be specific to the NLP */
    virtual bool ProcessOptions(const OptionsList& options,
                                const std::string& prefix)
    {
      return nlp_->ProcessOptions(options, prefix);
    }

    /** Method for creating the derived vector / matrix types.  The
     *  Hess_lagrangian_space pointer can be NULL if a quasi-Newton
     *  options is chosen. */
    virtual bool GetSpaces(SmartPtr<const VectorSpace>& x_space,
                           SmartPtr<const VectorSpace>& c_space,
                           SmartPtr<const VectorSpace>& d_space,
                           SmartPtr<const VectorSpace>& x_l_space,
                           SmartPtr<const MatrixSpace>& px_l_space,
                           SmartPtr<const VectorSpace>& x_u_space,
                           SmartPtr<const MatrixSpace>& px_u_space,
                           SmartPtr<const VectorSpace>& d_l_space,
                           SmartPtr<const MatrixSpace>& pd_l_space,
                           SmartPtr<const VectorSpace>& d_u_space,
                           SmartPtr<const MatrixSpace>& pd_u_space,
                           SmartPtr<const MatrixSpace>& Jac_c_space,
                           SmartPtr<const MatrixSpace>& Jac_d_space,
                           SmartPtr<const SymMatrixSpace>& Hess_lagrangian_space);

    /** Method for obtaining the bounds information */
    virtual bool GetBoundsInformation(const Matrix& Px_L,
                                      Vector& x_L,
                                      const Matrix& Px_U,
                                      Vector& x_U,
                                      const Matrix& Pd_L,
                                      Vector& d_L,
                                      const Matrix& Pd_U,
                                      Vector& d_U);

    /** Method for obtaining the starting point for all the
     *  iterates. ToDo it might not make sense to ask for initial
     *  values for v_L and v_U? */
    virtual bool GetStartingPoint(SmartPtr<Vector> x,
                                  bool need_x,
                                  SmartPtr<Vector> y_c,
                                  bool need_y_c,
                                  SmartPtr<Vector> y_d,
                                  bool need_y_d,
                                  SmartPtr<Vector> z_L,
                                  bool need_z_L,
                                  SmartPtr<Vector> z_U,
                                  bool need_z_U);

    /** Method for obtaining an entire iterate as a warmstart point.
     *  The incoming IteratesVector has to be filled.  This has not
     *  yet been implemented for this adapter. */
    virtual bool GetWarmStartIterate(IteratesVector& warm_start_iterate)
    {
      return nlp_->GetWarmStartIterate(warm_start_iterate);
    }
    //@}

    /** @name NLP evaluation routines (overload
     *  in derived classes. */
    //@{
    virtual bool Eval_f(const Vector& x, Number& f)
    {
      return nlp_->Eval_f(x, f);
    }

    virtual bool Eval_grad_f(const Vector& x, Vector& g_f)
    {
      return nlp_->Eval_grad_f(x, g_f);
    }

    virtual bool Eval_c(const Vector& x, Vector& c)
    {
      return nlp_->Eval_c(x, c);
    }

    virtual bool Eval_jac_c(const Vector& x, Matrix& jac_c)
    {
      return nlp_->Eval_jac_c(x, jac_c);
    }

    virtual bool Eval_d(const Vector& x, Vector& d);

    virtual bool Eval_jac_d(const Vector& x, Matrix& jac_d);

    virtual bool Eval_h(const Vector& x,
                        Number obj_factor,
                        const Vector& yc,
                        const Vector& yd,
                        SymMatrix& h);
    //@}

    /** @name NLP solution routines. Have default dummy
     *  implementations that can be overloaded. */
    //@{
    /** This method is called at the very end of the optimization.  It
     *  provides the final iterate to the user, so that it can be
     *  stored as the solution.  The status flag indicates the outcome
     *  of the optimization, where SolverReturn is defined in
     *  IpAlgTypes.hpp.  */
    virtual void FinalizeSolution(SolverReturn status,
                                  const Vector& x, const Vector& z_L,
                                  const Vector& z_U,
                                  const Vector& c, const Vector& d,
                                  const Vector& y_c, const Vector& y_d,
                                  Number obj_value,
                                  const IpoptData* ip_data,
                                  IpoptCalculatedQuantities* ip_cq);

    /** This method is called once per iteration, after the iteration
     *  summary output has been printed.  It provides the current
     *  information to the user to do with it anything she wants.  It
     *  also allows the user to ask for a premature termination of the
     *  optimization by returning false, in which case Ipopt will
     *  terminate with a corresponding return status.  The basic
     *  information provided in the argument list has the quantities
     *  values printed in the iteration summary line.  If more
     *  information is required, a user can obtain it from the IpData
     *  and IpCalculatedQuantities objects.  However, note that the
     *  provided quantities are all for the problem that Ipopt sees,
     *  i.e., the quantities might be scaled, fixed variables might be
     *  sorted out, etc.  The status indicates things like whether the
     *  algorithm is in the restoration phase...  In the restoration
     *  phase, the dual variables are probably not not changing. */
    virtual bool IntermediateCallBack(AlgorithmMode mode,
                                      Index iter, Number obj_value,
                                      Number inf_pr, Number inf_du,
                                      Number mu, Number d_norm,
                                      Number regularization_size,
                                      Number alpha_du, Number alpha_pr,
                                      Index ls_trials,
                                      const IpoptData* ip_data,
                                      IpoptCalculatedQuantities* ip_cq)
    {
      return nlp_->IntermediateCallBack(mode,iter, obj_value, inf_pr, inf_du,
                                        mu, d_norm, regularization_size,
                                        alpha_du, alpha_pr, ls_trials,
                                        ip_data, ip_cq);
    }
    //@}

    /** Routines to get the scaling parameters. These do not need to
     *  be overloaded unless the options are set for User scaling
     */
    //@{
    virtual void GetScalingParameters(
      const SmartPtr<const VectorSpace> x_space,
      const SmartPtr<const VectorSpace> c_space,
      const SmartPtr<const VectorSpace> d_space,
      Number& obj_scaling,
      SmartPtr<Vector>& x_scaling,
      SmartPtr<Vector>& c_scaling,
      SmartPtr<Vector>& d_scaling) const;
    //@}

    /** Method for obtaining the subspace in which the limited-memory
     *  Hessian approximation should be done.  This is only called if
     *  the limited-memory Hessian approximation is chosen.  Since the
     *  Hessian is zero in the space of all variables that appear in
     *  the problem functions only linearly, this allows the user to
     *  provide a VectorSpace for all nonlinear variables, and an
     *  ExpansionMatrix to lift from this VectorSpace to the
     *  VectorSpace of the primal variables x.  If the returned values
     *  are NULL, it is assumed that the Hessian is to be approximated
     *  in the space of all x variables.  The default instantiation of
     *  this method returns NULL, and a user only has to overwrite
     *  this method if the approximation is to be done only in a
     *  subspace. */
    virtual void
    GetQuasiNewtonApproximationSpaces(SmartPtr<VectorSpace>& approx_space,
                                      SmartPtr<Matrix>& P_approx)
    {
      nlp_->GetQuasiNewtonApproximationSpaces(approx_space, P_approx);
    }

    /** Accessor method to the original NLP */
    SmartPtr<NLP> nlp()
    {
      return nlp_;
    }

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default Constructor */
    NLPBoundsRemover();
    /** Copy Constructor */
    NLPBoundsRemover(const NLPBoundsRemover&);

    /** Overloaded Equals Operator */
    void operator=(const NLPBoundsRemover&);
    //@}

    /** Pointer to the original NLP */
    SmartPtr<NLP> nlp_;

    /** Pointer to the expansion matrix for the lower x bounds */
    SmartPtr<const Matrix> Px_l_orig_;

    /** Pointer to the expansion matrix for the upper x bounds */
    SmartPtr<const Matrix> Px_u_orig_;

    /** Pointer to the original d space */
    SmartPtr<const VectorSpace> d_space_orig_;

    /** Flag indicating whether twosided inequality constraints are
     *  allowed */
    bool allow_twosided_inequalities_;
  };

} // namespace Ipopt

#endif
