// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter, Frank E. Curtis         IBM    2009-06-12
//               (based on IpMc19TSymScalingMethod.hpp rev 699)

#ifndef __IPINEXACTTSYMSCALINGMETHOD_HPP__
#define __IPINEXACTTSYMSCALINGMETHOD_HPP__

#include "IpUtils.hpp"
#include "IpTSymScalingMethod.hpp"
#include "IpInexactCq.hpp"

namespace Ipopt
{

  /** Class for the method for computing scaling factors for symmetric
   *  matrices in triplet format, specifically for the inexaxct algorithm.
   *  The scaling is only considering the current slacks.
   */
  class InexactTSymScalingMethod: public TSymScalingMethod
  {
  public:
    /** @name Constructor/Destructor */
    //@{
    InexactTSymScalingMethod()
    {}

    virtual ~InexactTSymScalingMethod()
    {}
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Method for computing the symmetric scaling factors, given the
     *  symmtric matrix in triplet (MA27) format. */
    virtual bool ComputeSymTScalingFactors(Index n,
                                           Index nnz,
                                           const ipfint* airn,
                                           const ipfint* ajcn,
                                           const double* a,
                                           double* scaling_factors);
  private:
    /**@name Default Compiler Generated Methods (Hidden to avoid
     * implicit creation/calling).  These methods are not implemented
     * and we do not want the compiler to implement them for us, so we
     * declare them private and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Copy Constructor */
    InexactTSymScalingMethod(const InexactTSymScalingMethod&);

    /** Overloaded Equals Operator */
    void operator=(const InexactTSymScalingMethod&);

    /** Method to easily access Inexact calculated quantities */
    InexactCq& InexCq()
    {
      InexactCq& inexact_cq =
        static_cast<InexactCq&>(IpCq().AdditionalCq());
      DBG_ASSERT(dynamic_cast<InexactCq*>(&IpCq().AdditionalCq()));
      return inexact_cq;
    }

  };


} // namespace Ipopt

#endif
