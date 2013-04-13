// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-03-17

#ifndef __IPMC19TSYMSCALINGMETHOD_HPP__
#define __IPMC19TSYMSCALINGMETHOD_HPP__

#include "IpUtils.hpp"
#include "IpTSymScalingMethod.hpp"

namespace Ipopt
{

  /** Class for the method for computing scaling factors for symmetric
   *  matrices in triplet format, using MC19.
   */
  class Mc19TSymScalingMethod: public TSymScalingMethod
  {
  public:
    /** @name Constructor/Destructor */
    //@{
    Mc19TSymScalingMethod()
    {}

    virtual ~Mc19TSymScalingMethod()
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
    Mc19TSymScalingMethod(const Mc19TSymScalingMethod&);

    /** Overloaded Equals Operator */
    void operator=(const Mc19TSymScalingMethod&);
  };


} // namespace Ipopt

#endif
