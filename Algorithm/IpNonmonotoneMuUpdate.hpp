// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPNONMONOTONEMUUPDATE_HPP__
#define __IPNONMONOTONEMUUPDATE_HPP__

#include "IpMuUpdate.hpp"
#include "IpLineSearch.hpp"
#include "IpMuOracle.hpp"
#include "IpFilter.hpp"

namespace Ipopt
{

  /** Non-monotone mu update.
   */
  class NonmonotoneMuUpdate : public MuUpdate
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor */
    NonmonotoneMuUpdate(const SmartPtr<LineSearch>& linesearch,
                        const SmartPtr<MuOracle>& free_mu_oracle,
                        const SmartPtr<MuOracle>& fix_mu_oracle=NULL);
    /** Default destructor */
    virtual ~NonmonotoneMuUpdate();
    //@}

    /** Initialize method - overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Method for determining the barrier parameter for the next iteration.
     *  When the optimality error for the current barrier parameter is less than
     *  a tolerance, the barrier parameter is reduced, and the Reset method of the
     *  LineSearch object linesearch is called.
     *  TODO: MORE DETAILS HERE */
    virtual void UpdateBarrierParameter();

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
    NonmonotoneMuUpdate();

    /** Copy Constructor */
    NonmonotoneMuUpdate(const NonmonotoneMuUpdate&);

    /** Overloaded Equals Operator */
    void operator=(const NonmonotoneMuUpdate&);
    //@}

    /** @name Algorithmic parameters */
    //@{
    Number mu_max_;
    Number mu_min_;
    Number tau_min_;
    Number mu_safeguard_exp_;
    Number mu_safeguard_factor_; //ToDo don't need that?
    Number fixed_mu_avrg_factor_;
    /** ToDo the following should be combined with MonotoneMuUpdate */
    Number kappa_epsilon_;
    Number kappa_mu_;
    Number theta_mu_;
    Index nonmonotone_kkt_norm_;
    Index nonmonotone_kkt_centrality_;
    Index nonmonotone_kkt_balancing_term_;
    /** Flag indicating which globalization strategy should be used. */
    Index adaptive_globalization_;
    /** Maximal margin in filter (for adaptive_globalization = 3) */
    Number filter_max_margin_;
    /** Factor for filter margin */
    Number filter_margin_fact_;
    //@}

    /** @name Strategy objects */
    //@{
    /** Pointer to strategy object that is to be used for computing a
     *  suggested value of the barrier parameter in the free mu mode.
     */
    SmartPtr<MuOracle> free_mu_oracle_;
    /** Line search object of the Ipopt algorithm.  */
    SmartPtr<LineSearch> linesearch_;
    /** Pointer to strategy object that is to be used for computing a
     *  suggested value for the fixed mu mode.  If NULL, the current
     *  average complementarity is used.
     */
    SmartPtr<MuOracle> fix_mu_oracle_;
    //@}

    /** Dual infeasibility at initial point.  A negative value means
     *  that this quantity has not yet been initialized. */
    Number init_dual_inf_;
    /** Primal infeasibility at initial point.  A negative value means
     *  that this quantity has not yet been initialized. */
    Number init_primal_inf_;

    /** @name Methods and data defining the outer globalization
     *  strategy (might be a strategy object later). */
    //@{
    void InitializeFixedMuGlobalization();
    /** Check whether the point in the "current" fields offers
     *  sufficient reduction in order to remain in or switch to the
     *  free mu mode. */
    bool CheckSufficientProgress();
    /** Include the current point in internal memory to as accepted
     *  point */
    void RememberCurrentPointAsAccepted();
    /** Compute the value of the fixed mu that should be used in a new
     *  fixed mu phase.  This method is called at the beginning of a
     *  new fixed mu phase. */
    Number NewFixedMu();
    /** Compute value for the fraction-to-the-boundary parameter given
     *  mu in the monotone phase */
    Number Compute_tau_monotone(Number mu);

    /** Method for computing the 1-norm of the primal dual system at
     *  the current point.  The individual components (dual
     *  infeasibility, primal infeasibility, complementarity) are
     *  scaled to each other. */
    Number curr_norm_pd_system();

    /** Method for computing a lower safeguard bound for the barrier
     *  parameter.  For now, this is related to primal and dual
     *  infeasibility. */
    Number lower_mu_safeguard();

    /** Computer the currently largest reference value. */
    Number max_ref_val();

    /** Computer the currently smallest reference value. */
    Number min_ref_val();

    /** Maximal number of reference values (algorithmic parameter) */
    Index num_refs_max_;
    /** Values of the currently stored reference values (norm of pd
     *  equations) */
    std::list<Number> refs_vals_;
    /** Factor requested to reduce the reference values */
    Number refs_red_fact_;
    /** Flag indicating whether the barrier parameter should never be
     *  fixed (no globalization) */
    bool mu_never_fix_;

    /** Alternatively, we might also want to use a filter */
    Filter filter_;
    /** Flag indicating whether the most recent accepted step should
     *  be restored, when switching to the fixed mode. */
    bool restore_accepted_iterate_;
    //@}

    /** Flag indicating whether the problem has any inequality constraints */
    bool no_bounds_;
    /** Flag indicating whether no_bounds_ has been initialized */
    bool check_if_no_bounds_;

    /** @name Most recent accepted point in free mode, from which
     *  fixed mode should be started.
     */
    //@{
    SmartPtr<const Vector> accepted_x_;
    SmartPtr<const Vector> accepted_s_;
    SmartPtr<const Vector> accepted_y_c_;
    SmartPtr<const Vector> accepted_y_d_;
    SmartPtr<const Vector> accepted_z_L_;
    SmartPtr<const Vector> accepted_z_U_;
    SmartPtr<const Vector> accepted_v_L_;
    SmartPtr<const Vector> accepted_v_U_;
    //@}

  };

} // namespace Ipopt

#endif
