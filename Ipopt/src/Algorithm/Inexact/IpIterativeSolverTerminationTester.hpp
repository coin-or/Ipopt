// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2008-09-19

#ifndef __IPITERATIVESOLVERTERMINATIONTESTER_HPP__
#define __IPITERATIVESOLVERTERMINATIONTESTER_HPP__

#include "IpAlgStrategy.hpp"
#include "IpInexactCq.hpp"

namespace Ipopt
{

  /** This base class is for the termination tests for the iterative
   *  linear solver in the inexact version of Ipopt.
   */
  class IterativeSolverTerminationTester: public AlgorithmStrategyObject
  {
  public:
    /** Enum to report result of termination test */
    enum ETerminationTest {
      /** The current solution is not yet good enough */
      CONTINUE,
      /** Termination Test 1 is satisfied */
      TEST_1_SATISFIED,
      /** Termination Test 2 is satisfied */
      TEST_2_SATISFIED,
      /** Termination Test 3 is satisfied */
      TEST_3_SATISFIED,
      /** Hessian matrix should be modified */
      MODIFY_HESSIAN,
      /** Some other termination criterion satisfied */
      OTHER_SATISFIED
    };

    /** @name /Destructor */
    //@{
    /** Default constructor
     */
    IterativeSolverTerminationTester()
    {}

    /** Default destructor */
    virtual ~IterativeSolverTerminationTester()
    {}
    //@}

    /* overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix) = 0;


    /** Method for initializing for the next iterative solve.  This
     *  must be call before the test methods are called. */
    virtual bool InitializeSolve() = 0;

    /** This method checks if the current soltion of the iterative
     *  linear solver is good enough (by returning the corresponding
     *  satisfied termination test), or if the Hessian should be
     *  modified.  The input is the dimension of the augmented system,
     *  the current solution vector of the augmented system, the
     *  current residual vector. */
    virtual ETerminationTest TestTermination(Index ndim, const Number* sol,
        const Number* resid, Index iter,
        Number norm2_rhs) = 0;

    /** This method can be called after the Solve is over and we can
     *  delete anything that has been allocated to free memory. */
    virtual void Clear() = 0;

    /** An easy way to get the journalist if accessed from the outside */
    const Journalist& GetJnlst() const
    {
      return Jnlst();
    }

    /** Return the number of iterative solver iteration from the most
     *  recent solve */
    virtual Index GetSolverIterations() const = 0;

  protected:
    /** Method for copying a long augmented system array into Vectors
     *  in Ipopt notation */
    void GetVectors(Index ndim, const Number* array,
                    SmartPtr<const Vector>& comp_x,
                    SmartPtr<const Vector>& comp_s,
                    SmartPtr<const Vector>& comp_c,
                    SmartPtr<const Vector>& comp_d);

    /** Method to easily access Inexact data */
    InexactData& InexData()
    {
      InexactData& inexact_data =
        static_cast<InexactData&>(IpData().AdditionalData());
      DBG_ASSERT(dynamic_cast<InexactData*>(&IpData().AdditionalData()));
      return inexact_data;
    }

    /** Method to easily access Inexact calculated quantities */
    InexactCq& InexCq()
    {
      InexactCq& inexact_cq =
        static_cast<InexactCq&>(IpCq().AdditionalCq());
      DBG_ASSERT(dynamic_cast<InexactCq*>(&IpCq().AdditionalCq()));
      return inexact_cq;
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
    IterativeSolverTerminationTester& operator=(const IterativeSolverTerminationTester&);
    //@}
  };

} // namespace Ipopt

#endif
