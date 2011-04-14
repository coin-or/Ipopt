// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-14


#ifndef __ASASBACKSOLVER_HPP__
#define __ASASBACKSOLVER_HPP__

#include "IpAlgStrategy.hpp"
#include "IpIteratesVector.hpp"

namespace Ipopt
{

  class SensBacksolver : public AlgorithmStrategyObject
  {

    /** This class is the interface to all backsolvers that may
     *  be used for the sIPOPT. */
  public:
    SensBacksolver()
    {
    }

    virtual ~SensBacksolver()
    {
    }

    virtual bool Solve(SmartPtr<IteratesVector> delta_lhs, SmartPtr<const IteratesVector> delta_rhs)=0;

  };

}

#endif
