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
    bool IsAcceptableToCurrentIterate(Number trial_barr, Number trial_theta) const;

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

    /** Method to calculate alpha_min (minimum alpha before going to restoration
     */
    Number CalculateAlphaMin();

    /** Augment the filter used on the current values of the barrier
     *  objective function and the contraint violation */
    void AugmentFilter();

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
                                  SmartPtr<const Vector>& actual_delta_x,
                                  SmartPtr<const Vector>& actual_delta_s,
                                  SmartPtr<const Vector>& actual_delta_y_c,
                                  SmartPtr<const Vector>& actual_delta_y_d,
                                  SmartPtr<const Vector>& actual_delta_z_L,
                                  SmartPtr<const Vector>& actual_delta_z_U,
                                  SmartPtr<const Vector>& actual_delta_v_L,
                                  SmartPtr<const Vector>& actual_delta_v_U);

    /** Try higher order corrector (for fast local convergence).  In
     *  contrast to a second order correction step, which tries to
     *  make an unacceptable point acceptable by improving constraint
     *  violation, this corrector step is tried even if the regular
     *  primal-dual step is acceptable.
     */
    bool TryCorrector(Number alpha_primal_test,
                      Number& alpha_primal,
                      SmartPtr<const Vector>& actual_delta_x,
                      SmartPtr<const Vector>& actual_delta_s,
                      SmartPtr<const Vector>& actual_delta_y_c,
                      SmartPtr<const Vector>& actual_delta_y_d,
                      SmartPtr<const Vector>& actual_delta_z_L,
                      SmartPtr<const Vector>& actual_delta_z_U,
                      SmartPtr<const Vector>& actual_delta_v_L,
                      SmartPtr<const Vector>& actual_delta_v_U);

    /** Perform magic steps.  Take the current values of the slacks in
     *  trial and replace them by better ones that lead to smaller
     *  values of the barrier function and less constraint
     *  violation. */
    void PerformMagicStep();

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
     *  the equality constraint multipliers. */
    bool dual_alpha_for_y_;

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

    SmartPtr<RestorationPhase> resto_phase_;
    SmartPtr<PDSystemSolver> pd_solver_;
  };

} // namespace Ipopt

#endif
