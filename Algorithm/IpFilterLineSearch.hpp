// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPFILTERLINESEARCH_HPP__
#define __IPFILTERLINESEARCH_HPP__

#include "IpUtils.hpp"
#include "IpFilter.hpp"
#include "IpLineSearch.hpp"
#include "IpRestoPhase.hpp"
#include "IpPDSystemSolver.hpp"

namespace Ipopt
{

  /** Filter line search.  This class implements the filter line
   *  search procedure. 
   */
  class FilterLineSearch : public LineSearch
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor.  The PDSystemSolver object only needs to be
     *  provided (i.e. not NULL) if second order correction is to be
     *  used. */
    FilterLineSearch(const SmartPtr<RestorationPhase>& resto_phase,
                     const SmartPtr<PDSystemSolver>& pd_solver
                    );

    /** Default destructor */
    virtual ~FilterLineSearch();
    //@}

    /** InitializeImpl - overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Perform the line search.  It is assumed that the search
     *  direction is computed in the data object.
     */
    virtual void FindAcceptableTrialPoint();

    /** Reset the line search.
     *  This function should be called if all previous information
     *  should be discarded when the line search is performed the
     *  next time.  For example, this method should be called if
     *  the barrier parameter is changed.
     */
    virtual void Reset();

    /** Set flag indicating whether a very rigorous line search should
     *  be performed.  If this flag is set to true, the line search
     *  algorithm might decide to abort the line search and not to
     *  accept a new iterate.  If the line search decided not to
     *  accept a new iterate, the return value of
     *  CheckSkippedLineSearch() is true at the next call.  For
     *  example, in the non-monotone barrier parameter update
     *  procedure, the filter algorithm should not switch to the
     *  restoration phase in the free mode; instead, the algorithm
     *  should swtich to the fixed mode.
     */
    virtual void SetRigorousLineSearch(bool rigorous)
    {
      rigorous_ = rigorous;
    }

    /** Check if the line search procedure didn't accept a new iterate
     *  during the last call of FindAcceptableTrialPoint().
     *  
     */
    virtual bool CheckSkippedLineSearch()
    {
      return skipped_line_search_;
    }

    /**@name Trial Point Accepting Methods. Used internally to check certain
     * acceptability criteria and used externally (by the restoration phase
     * convergence check object, for instance)
     */
    //@{
    /** Checks if a trial point is acceptable to the current iterate */
    bool IsAcceptableToCurrentIterate(Number trial_barr, Number trial_theta,
                                      bool called_from_restoration=false) const;

    /** Checks if a trial point is acceptable to the current filter */
    bool IsAcceptableToCurrentFilter(Number trial_barr, Number trial_theta) const;
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
    FilterLineSearch(const FilterLineSearch&);

    /** Overloaded Equals Operator */
    void operator=(const FilterLineSearch&);
    //@}

    /** @name Filter information */
    //@{
    /** Upper bound on infeasibility */
    Number theta_max_;
    Number theta_max_fact_;

    /** Infeasibility switching bound */
    Number theta_min_;
    Number theta_min_fact_;

    //@}

    /** Method for checking if the current step size satisfies the
     *  f-type switching condition.  Here, we use the search direction
     *  stored in ip_data
     */
    bool IsFtype(Number alpha_primal_test);

    /** Method for checking the Armijo condition, given a trial step
     *  size.  The test uses the search direction stored in ip_data,
     *  and the values of the functions at the trial point in ip_data.
     */
    bool ArmijoHolds(Number alpha_primal_test);

    /** Method to calculate alpha_min (minimum alpha before going to
     *  restoration
     */
    Number CalculateAlphaMin();

    /** Augment the filter used on the current values of the barrier
     *  objective function and the contraint violation */
    void AugmentFilter();

    /** Method performing the backtracking line search.  The return
     *  value indicates if the step acceptance criteria are met.  If
     *  the watchdog is active, only one trial step is performed (and
     *  the trial values are set accordingly). */
    bool DoBacktrackingLineSearch(bool skip_first_trial_point,
                                  Number& alpha_primal,
                                  bool& corr_taken,
                                  bool& soc_taken,
                                  Index& n_steps,
                                  bool& evaluation_error,
				  SmartPtr<IteratesVector>& actual_delta);

    /** Method for starting the watch dog.  Set all appropriate fields
     *  accordingly */
    void StartWatchDog();

    /** Method for stopping the watch dog.  Set all appropriate fields
     *  accordingly. */
    void StopWatchDog(SmartPtr<IteratesVector>& actual_delta);

    /** Method for checking if current trial point is acceptable.
     *  It is assumed that the delta information in ip_data is the
     *  search direction used in criteria.  The primal trial point has
     *  to be set before the call.
     */
    bool CheckAcceptabilityOfTrialPoint(Number alpha_primal);

    /** Check comparison "lhs <= rhs", using machine precision based on BasVal */
    //ToDo This should probably not be a static member function if we want to
    //     allow for different relaxation parameters values
    static bool Compare_le(Number lhs, Number rhs, Number BasVal);

    /** Method for setting the dual variables in the trial fields in
     *  IpData, given the search direction.  The step size for the
     *  bound multipliers is alpha_dual (the fraction-to-the-boundary
     *  step size), and the step size for the equality constraint 
     *  multipliers depends on the choice of alpha_for_y. */
    void PerformDualStep(Number alpha_primal,
			 Number alpha_dual,
			 SmartPtr<IteratesVector>& delta);

    /** Try a step for the soft restoration phase and check if it is
     *  acceptable.  The step size is identical for all variables.  A
     *  point is accepted if it is acceptable for the original filter
     *  (in which case satisfies_original_filter = true on return), or
     *  if the primal-dual system error was decrease by at least the
     *  factor resto_pderror_reduction_factor_.  The return value is
     *  true, if the trial point was acceptable. */
    bool TrySoftRestoStep(SmartPtr<IteratesVector>& actual_delta,
                          bool &satisfies_original_filter);

    /** Try a second order correction for the constraints.  If the
     *  first trial step (with incoming alpha_primal) has been reject,
     *  this tries up to max_soc_ second order corrections for the
     *  constraints.  Here, alpha_primal_test is the step size that
     *  has to be used in the filter acceptance tests.  On output
     *  actual_delta_... has been set to the steps including the
     *  second order correction if it has been accepted, otherwise it
     *  is unchanged.  If the SOC step has been accepted, alpha_primal
     *  has the fraction-to-the-boundary value for the SOC step on output.
     *  The return value is true, if an SOC step has been accepted.
     */
    bool TrySecondOrderCorrection(Number alpha_primal_test,
                                  Number& alpha_primal,
				  SmartPtr<IteratesVector>& actual_delta);

    /** Try higher order corrector (for fast local convergence).  In
     *  contrast to a second order correction step, which tries to
     *  make an unacceptable point acceptable by improving constraint
     *  violation, this corrector step is tried even if the regular
     *  primal-dual step is acceptable.
     */
    bool TryCorrector(Number alpha_primal_test,
                      Number& alpha_primal,
		      SmartPtr<IteratesVector>& actual_delta);

    /** Perform magic steps.  Take the current values of the slacks in
     *  trial and replace them by better ones that lead to smaller
     *  values of the barrier function and less constraint
     *  violation. */
    void PerformMagicStep();

    /** Detect if the search direction is too small.  This should be
     *  true if the search direction is so small that if makes
     *  numerically no difference. */
    bool DetectTinyStep();

    /** Store current iterate as acceptable point */
    void StoreAcceptablePoint();

    /** Restore acceptable point into the current fields of IpData if
     *  found. Returns true if such as point is available. */
    bool RestoreAcceptablePoint();

    /** @name Parameters for the filter algorithm.  Names as in the paper */
    //@{
    /** \f$ \eta_{\varphi} \f$ */
    Number eta_phi_;
    /** \f$ \delta \f$ */
    Number delta_;
    /** \f$ s_{\varphi} \f$ */
    Number s_phi_;
    /** \f$ s_{\Theta} \f$ */
    Number s_theta_;
    /** \f$ \gamma_{\varphi} \f$ */
    Number gamma_phi_;
    /** \f$ \gamma_{\Theta} \f$ */
    Number gamma_theta_;
    /** \f$ \gamma_{\alpha} \f$ */
    Number alpha_min_frac_;
    /** factor by which search direction is to be shortened if trial
     *  point is rejected. */
    Number alpha_red_factor_;
    /** Maximal number of second order correction steps */
    Index max_soc_;
    /** Required reduction in constraint violation before trying
     *  multiple second order correction steps \f$ \kappa_{soc}\f$.
     */
    Number kappa_soc_;
    /** Maximal increase in objective function in orders of magnitute
     *  (log10).  If the log10(barrier objective function) is
     *  increased more than this compared to the current point, the
     *  trial point is rejected. */
    Number obj_max_inc_;
    /** Flag indicating whether the dual step size is to be used for
     *  the equality constraint multipliers. If 0, the primal step
     *  size is used, if 1 the dual step size, and if 2, the minimum
     *  of both. */
    Index alpha_for_y_;

    /** Flag indicating whether magic steps should be used. */
    bool magic_steps_;
    /** Type of corrector steps that should be tried. */
    Index corrector_type_;
    /** Flag indicating whether the line search should always accept
     *  the full (fraction-to-the-boundary) step. */
    bool ls_always_accept_;
    /** parameter in heurstic that determines whether corrector step
    should be tried. */
    Number corrector_compl_avrg_red_fact_;
    /** Flag indicating whether the corrector should be skipped in an
     *  iteration in which negative curvature is detected */
    bool skip_corr_if_neg_curv_;
    /** Flag indicating whether the corrector should be skipped in the
     *  fixed mu mode. */
    bool skip_corr_if_fixed_mode_;
    /** Indicates whether problem can be expected to be infeasible.
     *  This will trigger requesting a tighter reduction in
     *  infeasibility the first time the restoration phase is
     *  called. */
    bool expect_infeasible_problem_;
    /** Tolerance on constraint violation for
     *  expect_infeasible_problem heuristic.  If the constraint
     *  violation becomes that than this value, the heuristic is
     *  disabled for the rest of the optimization run. */
    Number expect_infeasible_problem_ctol_;
    /** Reduction factor for the restoration phase that accepts steps
     *  reducing the optimality error ("soft restoration phase"). If
     *  0., then this restoration phase is not enabled. */
    Number resto_pderror_reduction_factor_;

    /** Tolerance for detecting tiny steps. */
    Number tiny_step_tol_;

    /** Number of watch dog trial steps. */
    Index watch_dog_trial_iter_max_;
    /** Number of shortened iterations that trigger the watchdog. */
    Index watch_dog_shortened_iter_trigger_;

    /** Acceptable tolerance for the problem to terminate earlier if
     *  algorithm seems stuck or cycling */
    Number acceptable_tol_;
    /** Maximum number of iterations with acceptable level of accuracy
     *  and full steps, after which the algorithm terminates.  If 0,
     *  this heuristic is disabled. */
    Index acceptable_iter_max_;
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
    /** Flag indicating if the watchdog is active */
    bool in_watch_dog_;
    /** Counter for shortened iterations. */
    Index watch_dog_shortened_iter_;
    /** Counter for watch dog iterations */
    Index watch_dog_trial_iter_;
    /** Step size for Armijo test in watch dog */
    Number watch_dog_alpha_primal_test_;
    /** Constraint violation at reference point */
    Number watch_dog_theta_;
    /** Barrier objective function at reference point */
    Number watch_dog_barr_;
    /** Barrier gradient transpose search direction at reference point */
    Number watch_dog_gradBarrTDelta_;
    /** Watchdog reference iterate */
    SmartPtr<const IteratesVector> watch_dog_iterate_;
    /** Watchdog search direction at reference point */
    SmartPtr<const IteratesVector> watch_dog_delta_;
    //@}

    /** @name Storage for last iterate that satisfies the acceptable
     *  level of optimality error. */
    //@{
    SmartPtr<const IteratesVector> acceptable_iterate_;
    //@}

    /** Filter with entries */
    Filter filter_;

    /** Flag indicating whether the line search is to be performed
     * robust (usually this is true, unless SetRigorousLineSearch is
     * called with false).
     */
    bool rigorous_;

    /** Flag indicating whether no acceptable trial point was found
     *  during last line search. */
    bool skipped_line_search_;

    /** Flag indicating whether we are currently in the "soft"
     *  restoration phase mode, in which steps are accepted if they
     *  reduce the optimality error (see
     *  resto_pderror_reduction_factor) */
    bool in_soft_resto_phase_;

    /** Counter for the number of successive iterations in which the
     *  full step was not accepted. */
    Index count_successive_shortened_steps_;

    /** Counter for the number of successive iterations in which the
     *  nlp error was below the acceptable tolerance and a full step
     *  was accepted. */
    Index count_acceptable_iter_;

    /** Flag indicating if a tiny step was detected in previous
     *  iteration */
    bool tiny_step_last_iteration_;

    SmartPtr<RestorationPhase> resto_phase_;
    SmartPtr<PDSystemSolver> pd_solver_;
  };

} // namespace Ipopt

#endif
