// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2008-09-19

#ifndef __IPINEXACTNORMALTERMINATIONTESTER_HPP__
#define __IPINEXACTNORMALTERMINATIONTESTER_HPP__

#include "IpIterativeSolverTerminationTester.hpp"

namespace Ipopt
{

  /** This class implements the termination tests for the primal-dual
   *  system.
   */
  class InexactNormalTerminationTester: public IterativeSolverTerminationTester
  {
  public:
    /** @name /Destructor */
    //@{
    /** Default constructor
     */
    InexactNormalTerminationTester();

    /** Default destructor */
    virtual ~InexactNormalTerminationTester();
    //@}

    /* overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Methods for IpoptType */
    //@{
    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);
    //@}

    /** Method for initializing for the next iterative solve.  This
     *  must be call before the test methods are called. */
    virtual bool InitializeSolve();

    /** This method checks if the current soltion of the iterative
     *  linear solver is good enough (by returning the corresponding
     *  satisfied termination test), or if the Hessian should be
     *  modified.  The input is the dimension of the augmented system,
     *  the current solution vector of the augmented system, the
     *  current residual vector. */
    virtual ETerminationTest TestTerminaion(Index ndim, const Number* sol,
                                            const Number* resid, Index iter,
                                            Number norm2_rhs);

    /** This method can be called after the Solve is over and we can
     *  delete anything that has been allocated to free memory. */
    virtual void Clear();

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Overloaded Equals Operator */
    InexactNormalTerminationTester& operator=(const InexactNormalTerminationTester&);
    //@}

  };

} // namespace Ipopt

#endif
