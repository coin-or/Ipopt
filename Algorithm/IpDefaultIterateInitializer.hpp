// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter              IBM    2004-09-24

#ifndef __IPDEFAULTITERATEINITIALIZER_HPP__
#define __IPDEFAULTITERATEINITIALIZER_HPP__

#include "IpUtils.hpp"
#include "IpIterateInitializer.hpp"
#include "IpEqMultCalculator.hpp"

namespace Ipopt
{

  /** Class implementing the default initialization procedure (based
   *  on user options) for the iterates.  It is used at the very
   *  beginning of the optimization for determine the starting point
   *  for all variables.
   */
  class DefaultIterateInitializer: public IterateInitializer
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor.  If eq_mult_calculator is not NULL, it will be
     *  used to compute the initial values for equality constraint
     *  multipliers. */
    DefaultIterateInitializer
    (const SmartPtr<EqMultiplierCalculator>& eq_mult_calculator);

    /** Default destructor */
    virtual ~DefaultIterateInitializer()
    {}
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Compute the initial iterates and set the into the curr field
     *  of the ip_data object. */
    virtual bool SetInitialIterates();

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
    DefaultIterateInitializer();

    /** Copy Constructor */
    DefaultIterateInitializer(const DefaultIterateInitializer&);

    /** Overloaded Equals Operator */
    void operator=(const DefaultIterateInitializer&);
    //@}

    /**@name Algorithmic Parameters */
    //@{
    /** Parameters for bumping x0 */
    Number bound_push_;
    /** Parameters for bumping x0 */
    Number bound_frac_;
    /** If max-norm of the initial equality constraint multiplier
     *  estimate is larger than this, the initial y_* variables are
     *  set to zero. */
    Number laminitmax_;
    /** Initial value for all bound mulitpliers. */
    Number boundmultinitval_;
    //@}

    /** object to be used for the initialization of the equality
     *  constraint multipliers. */
    SmartPtr<EqMultiplierCalculator> eq_mult_calculator_;
  };

} // namespace Ipopt

#endif
