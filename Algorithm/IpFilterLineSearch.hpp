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
    /** Class for one filter entry. */
    class FilterEntry
    {
    public:
      /**@name Constructors/Destructors */
      //@{
      /** Constructor with the two components and the current iteration count */
      FilterEntry(Number phi, Number theta, Index iter);

      /** Default Destructor */
      ~FilterEntry();
      //@}

      /** Check acceptability of pair (phi,theta) with respect
       *  to this filter entry.  Returns true, if pair is acceptable.
       */
      bool Acceptable(Number phi, Number theta) const
      {
        return Compare_le(phi, phi_, phi_) ||
               Compare_le(theta, theta_, theta_);
      }

      /** Check if this entry is dominated by pair (phi,theta).
       * Returns true, if this entry is dominated.
       */
      bool Dominated(Number phi, Number theta) const
      {
        return (phi<=phi_ && theta<=theta_);
      }

      /** @name Accessor functions */
      //@{
      Number phi() const
      {
        return phi_;
      }
      Number theta() const
      {
        return theta_;
      }
      Index iter() const
      {
        return iter_;
      }
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
      /** Default Constructor */
      FilterEntry();
      /** Copy Constructor */
      FilterEntry(const FilterEntry&);

      /** Overloaded Equals Operator */
      void operator=(const FilterEntry&);
      //@}

      /** entry value for objective function */
      const Number phi_;
      /** entry value for constraint violation */
      const Number theta_;
      /** iteration number in which this entry was added to filter */
      const Index iter_;
    };

    /** Class for the filter.  This class contains all filter entries.
     *  The entries are stored as the corner point, including the
     *  margin. */
    class Filter
    {
    public:
      /**@name Constructors/Destructors */
      //@{
      /** Default Constructor */
      Filter()
      {}
      /** Default Destructor */
      ~Filter()
      {
        Clear();
      }
      //@}

      /** Check acceptability of pair (phi,theta) with respect
       *  to the filter.  Returns true, if pair is acceptable
       */
      bool Acceptable(Number phi, Number theta) const;

      /** Add filter entry for pair (phi,theta).  This will also
       *  delete all dominated entries in the current filter. */
      void AddEntry(Number phi, Number theta, Index iteration);

      /** Delete all filter entries */
      void Clear();

      /** Print current filter entries */
      void Print(const Journalist& jnlst);

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
      Filter(const Filter&);

      /** Overloaded Equals Operator */
      void operator=(const Filter&);
      //@}

      mutable std::list<FilterEntry*> filter_list_;
    };

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

    /** Resest the line search.
     *  This function should be called if all previous information
     *  should be discarded when the line search is performed the
     *  next time.  For example, this method should be called if
     *  the barrier parameter is changed.
     */
    virtual void Reset();

    /**@name Trial Point Accepting Methods. Used internally to check certain
     * acceptability criteria and used externally (by the restoration phase
     * convergence check object, for instance
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
    /** Number of entries in the filter */
    Index filter_size_;

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

    /** Compute search direction including a second order correction.
     *  Here, c_soc and d_minus_s_soc are assumed to be computed
     *  accoring to Eqn. (27) in paper.
     */
    void ComputeSecondOrderSearchDirection(const Vector& c_soc,
                                           const Vector& d_minus_s_soc,
                                           Vector& delta_soc_x,
                                           Vector& delta_soc_s,
                                           Vector& delta_soc_y_c,
                                           Vector& delta_soc_y_d,
                                           Vector& delta_soc_z_L,
                                           Vector& delta_soc_z_U,
                                           Vector& delta_soc_v_L,
                                           Vector& delta_soc_v_U);

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

    /** Flag indication whether magic steps should be used. */
    bool magic_steps_;
    //@}

    /** Filter with entries */
    Filter filter_;

    SmartPtr<RestorationPhase> resto_phase_;
    SmartPtr<PDSystemSolver> pd_solver_;
  };

} // namespace Ipopt

#endif
