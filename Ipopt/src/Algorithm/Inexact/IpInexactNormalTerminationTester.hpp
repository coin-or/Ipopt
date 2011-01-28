// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
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
    virtual ETerminationTest TestTermination(Index ndim, const Number* sol,
        const Number* resid, Index iter,
        Number norm2_rhs);

    /** This method can be called after the Solve is over and we can
     *  delete anything that has been allocated to free memory. */
    virtual void Clear();


    /** Return the number of iterative solver iteration from the most
     *  recent solve */
    virtual Index GetSolverIterations() const
    {
      return last_iter_;
    }

    /** Method for setting the normal problem objective function value
     *  at the Cauchy step.  This must be called by the Dogleg
     *  object. */
    void Set_c_Avc_norm_cauchy(Number c_Avc_norm_cauchy)
    {
      c_Avc_norm_cauchy_ = c_Avc_norm_cauchy;
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
    /** Overloaded Equals Operator */
    InexactNormalTerminationTester& operator=(const InexactNormalTerminationTester&);
    //@}

    /** @name Algorithmic options */
    //@{
    /** Desired reduction of residual */
    Number inexact_normal_tol_;
    /** Maximal number of iterative solve iterations */
    Index inexact_normal_max_iter_;
    /** Is set to true if the linear system is scaled via slacks. */
    bool requires_scaling_;
    //@}

    /** Value of normal problem objective function achived by the
     *  Cauchy step.  This must be set by the Dogleg step object. */
    Number c_Avc_norm_cauchy_;

    /** Last iterative solver iteration counter */
    Index last_iter_;
  };

} // namespace Ipopt

#endif
