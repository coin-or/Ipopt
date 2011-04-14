// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-06

#ifndef __SENSALGORITHM_HPP__
#define __SENSALGORITHM_HPP__

#include "IpAlgStrategy.hpp"
#include "SensStepCalc.hpp"
#include "SensMeasurement.hpp"
#include "SensSchurDriver.hpp"
#include "SensUtils.hpp"

namespace Ipopt
{

  class SensAlgorithm : public AlgorithmStrategyObject
  {
    /** This is the interface for the actual controller. It handles
     *  Data input to the controller (measurement) and returns controls */

  public:

    SensAlgorithm(std::vector< SmartPtr<SchurDriver> >& driver_vec,
		  SmartPtr<SensitivityStepCalculator> sens_step_calc,
		  SmartPtr<Measurement> measurement,
		  Index n_sens_steps);

    virtual ~SensAlgorithm();

    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Main loop: Wait for new measurement, Get new step, maybe deal with
     *  bounds,  see to it that everything happens in the required
     *  timeframe. */
    SensAlgorithmExitStatus Run();

  private:

    std::vector< SmartPtr<SchurDriver> > driver_vec_;
    SmartPtr<SensitivityStepCalculator> sens_step_calc_;
    SmartPtr<Measurement> measurement_;
    Index n_sens_steps_; // I think it is useful to state this number explicitly in the constructor and here.

  };
}

#endif
