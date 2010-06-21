// Copyright 2009, 2010 Hans Pirnay
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Date   : 2009-08-01

#ifndef __ASREDUCEDHESSIANCALCULATOR_HPP__
#define __ASREDUCEDHESSIANCALCULATOR_HPP__

#include "IpAlgStrategy.hpp"
#include "AsSchurData.hpp"
#include "AsPCalculator.hpp"

namespace Ipopt
{

  class ReducedHessianCalculator : public AlgorithmStrategyObject
  {
    /** This is the interface for the actual controller. It handles 
     *  Data input to the controller (measurement) and returns controls */
  public:
    ReducedHessianCalculator(SmartPtr<SchurData> hess_data,
			     SmartPtr<PCalculator> pcalc);

    virtual ~ReducedHessianCalculator();

    virtual bool InitializeImpl(const OptionsList& options,
				const std::string& prefix);

    /* This function computes the unscaled reduced hessian matrix */
    virtual bool ComputeReducedHessian();

  private:

    SmartPtr<SchurData> hess_data_;
    SmartPtr<PCalculator> pcalc_;
  };

}

#endif
