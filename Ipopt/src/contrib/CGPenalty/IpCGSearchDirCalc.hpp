// Copyright (C) 2005, 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2005-10-13
//
//           Lifeng Chen/Zaiwen Wen      Columbia Univ

#ifndef __IPCGSEARCHDIRCALC_HPP__
#define __IPCGSEARCHDIRCALC_HPP__

#include "IpSearchDirCalculator.hpp"
#include "IpPDSystemSolver.hpp"
#include "IpCGPenaltyCq.hpp"

namespace Ipopt
{

  /** Implementation of the search direction calculator that computes
   *  the Chen-Goldfarb step for the current barrier and penalty
   *  parameter.
   */
  class CGSearchDirCalculator : public SearchDirectionCalculator
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor */
    CGSearchDirCalculator(const SmartPtr<PDSystemSolver>& pd_solver);

    /** Default destructor */
    virtual ~CGSearchDirCalculator();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Method for computing the search direction.  If the penalty
     *  paraemeter has not yet been initialized, it is initialized
     *  now. The computed direction is stored in IpData().delta(). */
    virtual bool ComputeSearchDirection();

    /** Methods for IpoptType */
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
    /** Default Constructor */
    CGSearchDirCalculator();

    /** Copy Constructor */
    CGSearchDirCalculator(const CGSearchDirCalculator&);

    /** Overloaded Equals Operator */
    void operator=(const CGSearchDirCalculator&);
    //@}

    /** Method to easily access CGPenalty data */
    CGPenaltyData& CGPenData()
    {
      CGPenaltyData& cg_pen_data =
        static_cast<CGPenaltyData&>(IpData().AdditionalData());
      DBG_ASSERT(dynamic_cast<CGPenaltyData*>(&IpData().AdditionalData()));
      return cg_pen_data;
    }

    /** Method to easily access CGPenalty calculated quantities */
    CGPenaltyCq& CGPenCq()
    {
      CGPenaltyCq& cg_pen_cq =
        static_cast<CGPenaltyCq&>(IpCq().AdditionalCq());
      DBG_ASSERT(dynamic_cast<CGPenaltyCq*>(&IpCq().AdditionalCq()));
      return cg_pen_cq;
    }

    /** @name Algorithmic parameters */
    //@{
    /** safeguard factor for bound multipliers.  If value >= 1, then
     *  the dual variables will never deviate from the primal estimate
     *  by more than the factors kappa_sigma and 1./kappa_sigma.
     */
    Number penalty_init_min_;
    /** Maximal value for initial penalty parameter */
    Number penalty_init_max_;
    /** Maximal value for penalty parameters */
    Number penalty_max_;



    /**  parameters used in computation of line search penalty parameter and
     *    KKT perturbation parameters  **/
    Number pen_des_fact_;

    /** Algorithm type */
    bool penalty_backward_;

    /** parameters used to check if the fast direction can be
     *  used as the line search direction **/
    Number kappa_x_dis_;
    Number kappa_y_dis_;
    Number vartheta_;
    Number delta_y_max_;
    Number fast_des_fact_;
    Number pen_init_fac_;

    /** Flag indicating whether the fast Chen-Goldfarb direction
     *  should never be used */
    bool never_use_fact_cgpen_direction_;

    /** Counter for how many times the pen para is updated nonmonotonically*/
    Index nonmonotone_pen_update_counter_;
    //@}

    /** @name Strategy objects */
    //@{
    SmartPtr<PDSystemSolver> pd_solver_;
    //@}
  };

} // namespace Ipopt

#endif
