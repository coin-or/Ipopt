// Copyright (C) 2005 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2005-10-13

#ifndef __IPPDSEARCHDIRCALC_HPP__
#define __IPPDSEARCHDIRCALC_HPP__

#include "IpSearchDirCalculator.hpp"
#include "IpPDSystemSolver.hpp"

namespace Ipopt
{

  /** Implementation of the search direction calculator that computes
   *  the pure primal dual step for the current barrier parameter.
   */
  class PDSearchDirCalculator : public SearchDirectionCalculator
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor */
    PDSearchDirCalculator(const SmartPtr<PDSystemSolver>& pd_solver);

    /** Default destructor */
    virtual ~PDSearchDirCalculator();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Method for computing the search direction.  The computed
     *  direction is stored in IpData().delta(). */
    virtual void ComputeSearchDirection();

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
    PDSearchDirCalculator();

    /** Copy Constructor */
    PDSearchDirCalculator(const PDSearchDirCalculator&);

    /** Overloaded Equals Operator */
    void operator=(const PDSearchDirCalculator&);
    //@}

    /** @name Strategy objects */
    //@{
    SmartPtr<PDSystemSolver> pd_solver_;
    //@}
  };

} // namespace Ipopt

#endif
