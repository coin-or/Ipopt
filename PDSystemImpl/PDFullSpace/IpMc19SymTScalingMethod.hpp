// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPMC19SCALINGMETHOD_HPP__
#define __IPMC19SCLAINGMETHOD_HPP__

#include "IpUtils.hpp"
#include "IpAlgStrategy.hpp"

namespace Ipopt
{

  /** Class for the method for computing scaling factors for symmetric
   *  matrices in triplet format, based on MC19.
   */
  class Mc19SymTScalingMethod: public AlgorithmStrategyObject
  {
  public:
    /** @name Constructor/Destructor */
    //@{
    Mc19SymTScalingMethod()
    {}

    ~Mc19SymTScalingMethod()
    {}
    //@}

    /** overloaded from AlgorithmStrategyObject */
    bool InitializeImpl(const OptionsList& options,
                        const std::string& prefix);

    /** Method for computing the symmetric scaling factors, given the
     *  symmtric matrix in triplet (MA27) format. */
    bool ComputeSymTScalingFactors(Index n,
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
    Mc19SymTScalingMethod(const Mc19SymTScalingMethod&);

    /** Overloaded Equals Operator */
    void operator=(const Mc19SymTScalingMethod&);
  };


} // namespace Ipopt

#endif
