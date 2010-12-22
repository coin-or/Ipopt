// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter                 IBM    2008-09-11
//               derived file from IpPenaltyLSAcceptor.hpp (rev 019)

#ifndef __IPINEXACTLSACCEPTOR_HPP__
#define __IPINEXACTLSACCEPTOR_HPP__

#include "IpBacktrackingLSAcceptor.hpp"
#include "IpInexactCq.hpp"

namespace Ipopt
{

  /** Penalty function line search for the inexact step algorithm
   *  version.
   */
  class InexactLSAcceptor : public BacktrackingLSAcceptor
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor.  The PDSystemSolver object only needs to be
     *  provided (i.e. not NULL) if second order correction or
     *  corrector steps are to be used. */
    InexactLSAcceptor();

    /** Default destructor */
    virtual ~InexactLSAcceptor();
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
     *  indicates if we are currently in an active watchdog procedure.
     *  Here is where the penalty parameter is updated. */
    virtual void InitThisLineSearch(bool in_watchdog);

    /** Method that is called before the restoration phase is called.
     *  For now, we just terminate if this is called. */
    virtual void PrepareRestoPhaseStart();

    /** Method returning the lower bound on the trial step sizes. */
    virtual Number CalculateAlphaMin();

    /** Method for checking if current trial point is acceptable.
     *  It is assumed that the delta information in ip_data is the
     *  search direction used in criteria.  The primal trial point has
     *  to be set before the call.
     */
    virtual bool CheckAcceptabilityOfTrialPoint(Number alpha_primal);

    /** Try a second order correction for the constraints.  For the
     *  inexact version, this always returns false because a second
     *  order step is too expensive.
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
     *  the internal data should be updates.  alpha_primal_test is the
     *  value of alpha that has been used for in the acceptence test
     *  ealier. */
    virtual char UpdateForNextIteration(Number alpha_primal_test);

    /** Method for setting internal data if the watchdog procedure is
     *  started. */
    virtual void StartWatchDog();

    /** Method for setting internal data if the watchdog procedure is
     *  stopped. */
    virtual void StopWatchDog();

    /**@name Trial Point Accepting Methods. Used internally to check certain
     * acceptability criteria and used externally (by the restoration phase
     * convergence check object, for instance)
     */
    //@{
    /** Checks if a trial point is acceptable to the current iterate */
    bool IsAcceptableToCurrentIterate(Number trial_barr, Number trial_theta,
                                      bool called_from_restoration=false) const;
    //@}

    /** Method for updating the equality constraint multipliers */
    virtual Number ComputeAlphaForY(Number alpha_primal,
                                    Number alpha_dual,
                                    SmartPtr<IteratesVector>& delta);

    /** Method returning true of ComputeAlphaForY is implemented for
      *  this acceptor */
    virtual bool HasComputeAlphaForY() const
    {
      return true;
    }

    /** Methods for OptionsList */
    //@{
    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);
    //@}

  protected:
    /** Method to easily access Inexact data */
    InexactData& InexData()
    {
      InexactData& inexact_data =
        static_cast<InexactData&>(IpData().AdditionalData());
      DBG_ASSERT(dynamic_cast<InexactData*>(&IpData().AdditionalData()));
      return inexact_data;
    }

    /** Method to easily access Inexact calculated quantities */
    InexactCq& InexCq()
    {
      InexactCq& inexact_cq =
        static_cast<InexactCq&>(IpCq().AdditionalCq());
      DBG_ASSERT(dynamic_cast<InexactCq*>(&IpCq().AdditionalCq()));
      return inexact_cq;
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
    /** Copy Constructor */
    InexactLSAcceptor(const InexactLSAcceptor&);

    /** Overloaded Equals Operator */
    void operator=(const InexactLSAcceptor&);
    //@}

    /** Compute predicted reduction for given step size */
    Number CalcPred(Number alpha);

    /** Method for resetting the slacks to be satisfying the slack
    equality constraints without increasing the barrier
    function */
    void ResetSlacks();

    /** @name Parameters for the penalty function algorithm. */
    //@{
    /** Initial value of penalty parameter */
    Number nu_init_;
    /** Initial value of lower penalty parameter */
    Number nu_low_init_;
    /** Factor in update rule for lower penalty parameter */
    Number nu_low_fact_;
    /** Incrememt for penalty parameter */
    Number nu_inc_;
    /** \f$ \eta_{\varphi} \f$ */
    Number eta_;
    /** \f$ \rho \f$ */
    Number rho_;
    /** theta factor in Tangential Component Condition */
    Number tcc_theta_;
    /** Lower feasiblity bound to skip penalty parameter update */
    Number nu_update_inf_skip_tol_;
    /** Flag indicating whether the Curtis/Nocedal flexible penalty
     *  function should be used */
    bool flexible_penalty_function_;
    //@}

    /** @name Information related to watchdog procedure */
    //@{
    /** Constraint violation at the point with respect to which
     *  progress is to be made */
    Number reference_theta_;
    /** Barrier objective function at the point with respect to which
     *  progress is to be made */
    Number reference_barr_;
    /** Reference predicted reduction.  If positive, then it is used
     *  in watch dog. */
    Number reference_pred_;
    /** Constraint violation at reference point */
    Number watchdog_theta_;
    /** Barrier objective function at reference point */
    Number watchdog_barr_;
    /** Predicted reduction to be compared with in watch dog. */
    Number watchdog_pred_;
    //@}

    /** @name Penalty parameter */
    //@{
    /** Current value of the penalty parameter */
    Number nu_;
    /** Value of penalty parameter at beginning of the iteration. */
    Number last_nu_;
    /** Current lower value of the penalty parameter */
    Number nu_low_;
    /** Value of lower penalty parameter at beginning of the iteration */
    Number last_nu_low_;
    /** Step size threshold for activating step decomposition */
    Number inexact_decomposition_activate_tol_;
    /** Step size threshold for inactivating step decomposition */
    Number inexact_decomposition_inactivate_tol_;
    //@}

    /** Flag indicating if this is a termination test 2 iteration in
     *  which we just update the multipliers and skip the line
     *  search */
    bool in_tt2_;

    /** When called from the restoration phase, this is the required
     *  predicted reduction */
    Number resto_pred_;

    /** Flag indicating if the step was accepted only because of the lower
     *  penalty parameter.  This is for output only. */
    bool accepted_by_low_only_;
  };

} // namespace Ipopt

#endif
