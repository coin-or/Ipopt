// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPNONMONOTONEMUUPDATE_HPP__
#define __IPNONMONOTONEMUUPDATE_HPP__

#include "IpMuUpdate.hpp"
#include "IpLineSearch.hpp"
#include "IpMuOracle.hpp"

namespace Ipopt
{

  /** Non-monotone mu update.
   */
  class NonmonotoneMuUpdate : public MuUpdate
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor */
    NonmonotoneMuUpdate(const SmartPtr<LineSearch>& linesearch,
                        const SmartPtr<MuOracle>& mu_oracle);
    /** Default destructor */
    virtual ~NonmonotoneMuUpdate();
    //@}

    /** Initialize method - overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Method for determining the barrier parameter for the next iteration.
     *  When the optimality error for the current barrier parameter is less than
     *  a tolerance, the barrier parameter is reduced, and the Reset method of the
     *  LineSearch object linesearch is called.
     *  TODO: MORE DETAILS HERE */
    virtual void UpdateBarrierParameter();

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
    NonmonotoneMuUpdate();

    /** Copy Constructor */
    NonmonotoneMuUpdate(const NonmonotoneMuUpdate&);

    /** Overloaded Equals Operator */
    void operator=(const NonmonotoneMuUpdate&);
    //@}

    /** @name Algorithmic parameters */
    //@{
    Number mu_max_;
    Number mu_min_;
    Number tau_min_;
    //@}

    /** Pointer to the class that is to be used for computing a
     *  suggested value of the barrier parameter
     */
    SmartPtr<MuOracle> mu_oracle_;
    SmartPtr<LineSearch> linesearch_;
  };

} // namespace Ipopt

#endif
