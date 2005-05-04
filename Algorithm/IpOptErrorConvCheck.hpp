// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPOPTERRORCONVCHECK_HPP__
#define __IPOPTERRORCONVCHECK_HPP__

#include "IpUtils.hpp"
#include "IpConvCheck.hpp"

namespace Ipopt
{

  /** Brief Class Description.
   *  Detailed Class Description.
   */
  class OptimalityErrorConvergenceCheck : public ConvergenceCheck
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default Constructor */
    OptimalityErrorConvergenceCheck();

    /** Default destructor */
    virtual ~OptimalityErrorConvergenceCheck();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Overloaded convergence check */
    virtual ConvergenceStatus CheckConvergence();

  private:
    /**@name Default Compiler Generated Methods (Hidden to avoid
     * implicit creation/calling).  These methods are not implemented
     * and we do not want the compiler to implement them for us, so we
     * declare them private and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Copy Constructor */
    OptimalityErrorConvergenceCheck(const OptimalityErrorConvergenceCheck&);

    /** Overloaded Equals Operator */
    void operator=(const OptimalityErrorConvergenceCheck&);
    //@}

    /** @name Algorithmic parameters */
    //@{
    /** Maximal number of iterations */
    Index max_iterations_;
    //@}
  };

} // namespace Ipopt

#endif
