// Copyright (C) 2005, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
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
    virtual bool ComputeSearchDirection();

    /** Methods for IpoptType */
    //@{
    static void RegisterOptions(const SmartPtr<RegisteredOptions>& roptions);
    //@}

    /** Method to return the pd_solver for additional processing */
    SmartPtr<PDSystemSolver> PDSolver()
    {
      return pd_solver_;
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

    /** @name Algorithmic parameters */
    //@{
    /** Flag indicating that we trust that the steps from the linear
     *  solver are very good and that we don't need any residual
     *  checks */
    bool fast_step_computation_;
    /** Flag indicating if we want to do Mehrotras's algorithm.  This
     *  means that a number of options are ignored, or have to be set
     *  (or are automatically set) to certain values. */
    bool mehrotra_algorithm_;
    //@}

  };

} // namespace Ipopt

#endif
