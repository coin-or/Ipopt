// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPRESTOPHASE_HPP__
#define __IPRESTOPHASE_HPP__

#include "IpUtils.hpp"
#include "IpAlgStrategy.hpp"
#include "IpIpoptNLP.hpp"
#include "IpIpoptData.hpp"
#include "IpIpoptCalculatedQuantities.hpp"

namespace Ipopt
{

  /** Base class for different restoration phases.  The restoration
   *  phase is part of the FilterLineSearch. */
  class RestorationPhase : public AlgorithmStrategyObject
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default Constructor */
    RestorationPhase()
    {}
    /** Default Destructor */
    virtual ~RestorationPhase()
    {}
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix) = 0;

    /** Method called to perform restoration for the filter line
     *  search method.  Derived classes should not overload this
     *  method, but should overload the Impl version. */
    virtual bool PerformRestoration() = 0;

  private:
    /**@name Default Compiler Generated Methods (Hidden to avoid
     * implicit creation/calling).  These methods are not implemented
     * and we do not want the compiler to implement them for us, so we
     * declare them private and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Copy Constructor */
    RestorationPhase(const RestorationPhase&);

    /** Overloaded Equals Operator */
    void operator=(const RestorationPhase&);
    //@}
  };

} // namespace Ipopt

#endif
