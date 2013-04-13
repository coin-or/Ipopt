// Copyright (C) 2008, 2011 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2008-09-19

#ifndef __IPINEXACTPDTERMINATIONTESTER_HPP__
#define __IPINEXACTPDTERMINATIONTESTER_HPP__

#include "IpIterativeSolverTerminationTester.hpp"

namespace Ipopt
{

  /** This class implements the termination tests for the primal-dual
   *  system.
   */
  class InexactPDTerminationTester: public IterativeSolverTerminationTester
  {
  public:
    /** @name /Destructor */
    //@{
    /** Default constructor
     */
    InexactPDTerminationTester();

    /** Default destructor */
    virtual ~InexactPDTerminationTester();
    //@}

    /* overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Methods for IpoptType */
    //@{
    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);
    //@}

    /** Method for initializing for the next iterative solve.  This
     *  must be call before the test methods are called. */
    virtual bool InitializeSolve();

    /** This method checks if the current soltion of the iterative
     *  linear solver is good enough (by returning the corresponding
     *  satisfied termination test), or if the Hessian should be
     *  modified.  The input is the dimension of the augmented system,
     *  the current solution vector of the augmented system, the
     *  current residual vector. */
    virtual ETerminationTest TestTermination(Index ndim, const Number* sol,
        const Number* resid, Index iter,
        Number norm2_rhs);

    /** This method can be called after the Solve is over and we can
     *  delete anything that has been allocated to free memory. */
    virtual void Clear();

    /** Return the number of iterative solver iteration from the most
     *  recent solve */
    virtual Index GetSolverIterations() const
    {
      return last_iter_;
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
    /** Overloaded Equals Operator */
    InexactPDTerminationTester& operator=(const InexactPDTerminationTester&);
    //@}

    /** @name Algorithmic options */
    //@{
    /** Psi factor in the tangential component condition */
    Number tcc_psi_;
    /** theta factor in the tangential component condition */
    Number tcc_theta_;
    /** mu exponent when multiplied to theta in the tangential
     *  component condition */
    Number tcc_theta_mu_exponent_;
    /** zeta factor in the tangential component condition */
    Number tcc_zeta_;
    /** kappa_1 factor in termination test 1 */
    Number tt_kappa1_;
    /** kappa_2 factor in termination test 2 */
    Number tt_kappa2_;
    /** eps_2 constant in termination test 2 */
    Number tt_eps2_;
    /** eps_3 constant in termination test 3 */
    Number tt_eps3_;
    /** rho constant from penalty parameter update.  This is called
     *  \f$\tau_{\pi}\f$ in MIPS paper */
    Number rho_;
    /** Desired reduction of residual */
    Number inexact_desired_pd_residual_;
    /** Number of iterations allowed for desired pd residual */
    Index inexact_desired_pd_residual_iter_;
    /** Is set to true if the linear system is scaled via slacks. */
    bool requires_scaling_;
    //@}

    /** @name Quantities that are identical for all tests and can be
     *  precomputed */
    //@{
    SmartPtr<const Vector> curr_Av_c_;
    SmartPtr<const Vector> curr_Av_d_;
    Number c_norm_;
    Number c_plus_Av_norm_;
    Number v_norm_scaled_;
    SmartPtr<const Vector> curr_grad_barrier_obj_x_;
    SmartPtr<const Vector> curr_grad_barrier_obj_s_; // in original space
    SmartPtr<const Matrix> curr_jac_c_;
    SmartPtr<const Matrix> curr_jac_d_;
    SmartPtr<const Vector> curr_scaling_slacks_;
    SmartPtr<Vector> curr_nabla_phi_plus_ATy_x_;
    SmartPtr<Vector> curr_nabla_phi_plus_ATy_s_; // in scaled space
    Number curr_Av_norm_;
    Number curr_tt1_norm_;
    Number curr_tt2_norm_;
    SmartPtr<const Vector> curr_Wv_x_;
    SmartPtr<const Vector> curr_Wv_s_; // in original space
    bool try_tt2_;
    //@}

    /** @name Quantities from previous iteration required in the
    tests */
    //@{
    Number last_Av_norm_;
    Number last_tt1_norm_;
    //@}

    /** Last iterative solver iteration counter */
    Index last_iter_;
  };

} // namespace Ipopt

#endif
