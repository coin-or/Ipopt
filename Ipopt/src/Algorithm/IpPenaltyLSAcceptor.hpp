// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter                 IBM    2008-04-04
//               derived file from IpFilterLSAcceptor.hpp

#ifndef __IPPENALTYLSACCEPTOR_HPP__
#define __IPPENALTYLSACCEPTOR_HPP__

#include "IpBacktrackingLSAcceptor.hpp"
#include "IpPDSystemSolver.hpp"

namespace Ipopt
{

  /** Penalty function line search.  This class implements the penalty
   *  function line search procedure as proposed by Waltz, Morales,
   *  Nocedal, Orban.
   */
  class PenaltyLSAcceptor : public BacktrackingLSAcceptor
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor.  The PDSystemSolver object only needs to be
     *  provided (i.e. not NULL) if second order correction or
     *  corrector steps are to be used. */
    PenaltyLSAcceptor(const SmartPtr<PDSystemSolver>& pd_solver);

    /** Default destructor */
    virtual ~PenaltyLSAcceptor();
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

    /** Method returning the lower bound on the trial step sizes. */
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
     *  has to be used in the penalty function acceptance tests.  On
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
    PenaltyLSAcceptor(const PenaltyLSAcceptor&);

    /** Overloaded Equals Operator */
    void operator=(const PenaltyLSAcceptor&);
    //@}

    /** Compute predicted reduction for given step size */
    Number CalcPred(Number alpha);

    /** @name Parameters for the penalty function line search
     *  algorithm.  Names as in the filter paper */
    //@{
    /** Initial value of penalty parameter */
    Number nu_init_;
    /** Incrememt for penalty parameter */
    Number nu_inc_;
    /** \f$ \eta_{\varphi} \f$ */
    Number eta_;
    /** \f$ \rho \f$ */
    Number rho_;
    /** Maximal number of second order correction steps */
    Index max_soc_;
    /** Required reduction in constraint violation before trying
     *  multiple second order correction steps \f$ \kappa_{soc}\f$.
     */
    Number kappa_soc_;
    //@}

    /** @name Information related to watchdog procedure */
    //@{
    /** Constraint violation at the point with respect to which
     *  progress is to be made */
    Number reference_theta_;
    /** Barrier objective function at the point with respect to which
     *  progress is to be made */
    Number reference_barr_;
    /** Barrier gradient transpose search direction at the point with
     *  respect to which progress is to be made */
    Number reference_gradBarrTDelta_;
    /** Two-sided product of search direction with complete Hessian */
    Number reference_dWd_;
    /** Product of Jacobian of equality constraint with x direction */
    SmartPtr<const Vector> reference_JacC_delta_;
    /** Product of Jacobian of (d-s) constraint with search direction */
    SmartPtr<const Vector> reference_JacD_delta_;
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
    //@}

    /** When called from the restoration phase, this is the required
     *  predicted reduction */
    Number resto_pred_;

    /** @name Strategy objective that are used */
    //@{
    SmartPtr<PDSystemSolver> pd_solver_;
    //@}
  };

} // namespace Ipopt

#endif
