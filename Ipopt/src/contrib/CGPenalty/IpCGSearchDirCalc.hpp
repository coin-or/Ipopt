// Copyright (C) 2005, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpCGSearchDirCalc.hpp 551 2005-10-27 00:31:28Z andreasw $
//
// Authors:  Andreas Waechter            IBM    2005-10-13

#ifndef __IPCGSEARCHDIRCALC_HPP__
#define __IPCGSEARCHDIRCALC_HPP__

#include "IpSearchDirCalculator.hpp"
#include "IpPDSystemSolver.hpp"

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

    /** @name Algorithmic parameters */
    //@{
    /** safeguard factor for bound multipliers.  If value >= 1, then
     *  the dual variables will never deviate from the primal estimate
     *  by more than the factors kappa_sigma and 1./kappa_sigma.
     */
    Number penalty_init_min_;
    /** Maximal value for initial penalty parameter. */
    Number penalty_init_max_;
    /** Flag indicating whether the fast Chen-Goldfarb direction
     *  should never be used */
    bool never_use_fact_cgpen_direction_;
    //@}

    /** @name Strategy objects */
    //@{
    SmartPtr<PDSystemSolver> pd_solver_;
    //@}
  };

} // namespace Ipopt

#endif
