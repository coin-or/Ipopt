// Copyright 2009, 2010 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-16

#ifndef __ASSTDSTEPCALC_HPP__
#define __ASSTDSTEPCALC_HPP__

#include "AsSensStepCalc.hpp"
#include <vector>


namespace Ipopt
{

  class StdStepCalculator : public SensitivityStepCalculator
  {
  public: 
    StdStepCalculator();

    virtual ~StdStepCalculator();

    virtual bool InitializeImpl(const OptionsList& options,
				const std::string& prefix);

    /** This is the main algorithmic function of this class; It calculates
     *  a step using its SchurDriver, checks bounds, and returns it */
    virtual bool Step(DenseVector& delta_u, IteratesVector& sol);

    bool BoundCheck(IteratesVector& sol,
		    std::vector<Index>& x_bound_violations_idx,
		    std::vector<Number>& x_bound_violations_du);

  private:
    Number bound_eps_;
  };
}

#endif
