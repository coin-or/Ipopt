// Copyright (C) 2005, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpCGPenaltyLSAcceptor.hpp 551 2005-10-27 00:31:28Z andreasw $
//
// Authors:  Andreas Waechter                 IBM    2005-10-13

#ifndef __IPCGPENALTYLSACCEPTOR_HPP__
#define __IPCGPENALTYLSACCEPTOR_HPP__

#include "IpBacktrackingLSAcceptor.hpp"
#include "IpPDSystemSolver.hpp"

namespace Ipopt
{

  /** Line search acceptor, based on the Chen-Goldfarb penalty
   *  function approach. */
  class CGPenaltyLSAcceptor : public BacktrackingLSAcceptor
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor.  The PDSystemSolver object only needs to be
     *  provided (i.e. not NULL) if second order correction or
     *  corrector steps are to be used. */
    CGPenaltyLSAcceptor(const SmartPtr<PDSystemSolver>& pd_solver);

    /** Default destructor */
    virtual ~CGPenaltyLSAcceptor();
    //@}

    /** InitializeImpl - overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Reset the acceptor.
     *  This function should be called if all previous information
     *  should be discarded when the line search is performed the
     *  next time.  For example, this method should be called if
     *  the barrier parameter is changed.
     */
    virtual void Reset();

    /** Initialization for the next line search.  The flag in_watchdog
     *  indicates if we are currently in an active watchdog
     *  procedure. */
    virtual void InitThisLineSearch(bool in_watchdog);

    /** Method that is called before the restoration phase is called.
     *  Here, we can set up things that are required in the
     *  termination test for the restoration phase. */
    virtual void PrepareRestoPhaseStart();

    /** Method returning the lower bound on the trial step sizes.  If
     *  the backtracking procedure encounters a trial step size below
     *  this value after the first trial set, it swtiches to the
     *  (soft) restoration phase. */
    virtual Number CalculateAlphaMin();

    /** Method for checking if current trial point is acceptable.
     *  It is assumed that the delta information in ip_data is the
     *  search direction used in criteria.  The primal trial point has
     *  to be set before the call.
     */
    virtual bool CheckAcceptabilityOfTrialPoint(Number alpha_primal);

    /** Try a second order correction for the constraints.  If the
     *  first trial step (with incoming alpha_primal) has been reject,
     *  this tries up to max_soc_ second order corrections for the
     *  constraints.  Here, alpha_primal_test is the step size that
     *  has to be used in the merit function acceptance tests.  On
     *  output actual_delta_ has been set to the step including the
     *  second order correction if it has been accepted, otherwise it
     *  is unchanged.  If the SOC step has been accepted, alpha_primal
     *  has the fraction-to-the-boundary value for the SOC step on
     *  output.  The return value is true, if a SOC step has been
     *  accepted.
     */
    virtual bool TrySecondOrderCorrection(Number alpha_primal_test,
                                          Number& alpha_primal,
                                          SmartPtr<IteratesVector>& actual_delta);

    /** Try higher order corrector (for fast local convergence).  In
     *  contrast to a second order correction step, which tries to
     *  make an unacceptable point acceptable by improving constraint
     *  violation, this corrector step is tried even if the regular
     *  primal-dual step is acceptable.
     */
    virtual bool TryCorrector(Number alpha_primal_test,
                              Number& alpha_primal,
                              SmartPtr<IteratesVector>& actual_delta);

    /** Method for ending the current line search.  When it is called,
     *  the internal data should be updates, e.g., the penalty
     *  parameter might be updated.  alpha_primal_test is the value of
     *  alpha that has been used for in the acceptence test ealier. */
    virtual char UpdateForNextIteration(Number alpha_primal_test);

    /** Method for setting internal data if the watchdog procedure is
     *  started. */
    virtual void StartWatchDog();

    /** Method for setting internal data if the watchdog procedure is
     *  stopped. */
    virtual void StopWatchDog();

    /** Methods for OptionsList */
    //@{
    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);
    //@}

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Copy Constructor */
    CGPenaltyLSAcceptor(const CGPenaltyLSAcceptor&);

    /** Overloaded Equals Operator */
    void operator=(const CGPenaltyLSAcceptor&);
    //@}

    /** Check comparison "lhs <= rhs", using machine precision based on BasVal */
    //ToDo This should probably not be a static member function if we want to
    //     allow for different relaxation parameters values
    static bool Compare_le(Number lhs, Number rhs, Number BasVal);

    /** @name Parameters for the penalty function algorithm. */
    //@{
    /** Relaxation factor in the Armijo condition for the penalty function */
    Number eta_penalty_;
    /** Tolerance for infeasibility part in penalty parameter update
     *  rule. */
    Number penalty_update_infeasibility_tol_;
    /** Minimal tolerance for step part in penalty parameter update
     *  rule. */
    Number eta_min_;
    /** Tolerance for cimplementarity part in penalty parameter update
     *  rule. */
    Number penalty_update_compl_tol_;
    Number chi_hat_;
    Number chi_tilde_;
    Number chi_cup_;
    Number gamma_hat_;
    Number gamma_tilde_;
    Number penalty_max_;
    Number epsilon_c_;
    /** Maximal number of second order correction steps */
    Index max_soc_;
    /** Required reduction in constraint violation before trying
     *  multiple second order correction steps \f$ \kappa_{soc}\f$.
     */
    Number kappa_soc_;
    //@}

    /** Counter for increases of penalty parameter. */
    Index counter_penalty_updates_;

    /** eta parameter */
    Number curr_eta_;

    /** @name Information related to watchdog procedure */
    //@{
    /** Penalty function at the point with respect to which
     *  progress is to be made */
    Number reference_penalty_function_;
    /** Directional derivative of penalty function at the point with
     *  respect to which progress is to be made */
    Number reference_direct_deriv_penalty_function_;
    /** Penalty function at the point with respect to which
     *  progress is to be made (at watchdog point) */
    Number watchdog_penalty_function_;
    /** Directional derivative of penalty function at the point with
     *  respect to which progress is to be made (at watchdog point) */
    Number watchdog_direct_deriv_penalty_function_;
    /** Backup for the Chen-Goldfarb search direction (needed in the
     *  update rule for the penalty parameter */
    SmartPtr<const IteratesVector> watchdog_delta_cgpen_;
    //@}

    /** @name Strategy objective that are used */
    //@{
    SmartPtr<PDSystemSolver> pd_solver_;
    //@}
  };

} // namespace Ipopt

#endif
