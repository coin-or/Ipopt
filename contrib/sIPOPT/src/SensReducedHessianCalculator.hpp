// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-08-01

#ifndef __ASREDUCEDHESSIANCALCULATOR_HPP__
#define __ASREDUCEDHESSIANCALCULATOR_HPP__

#include "IpAlgStrategy.hpp"
#include "SensSchurData.hpp"
#include "SensPCalculator.hpp"

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

    /** Pointer to Schurdata object holding the indices for selecting the free variables */
    SmartPtr<SchurData> hess_data_;

    /** Pointer to the P Calculator object that returns the reduced hessian matrix */
    SmartPtr<PCalculator> pcalc_;

    /** True, if option rh_eigendecomp was set to yes */
    bool compute_eigenvalues_;
  };

}

#endif
