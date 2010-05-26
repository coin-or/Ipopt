// Copyright 2009 Hans Pirnay
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Date   : 2009-05-06

#ifndef __ASASNMPCONTROLLER_HPP__
#define __ASASNMPCONTROLLER_HPP__

#include "IpAlgStrategy.hpp"
#include "AsSensStepCalc.hpp"
#include "AsMeasurement.hpp"
#include "AsSchurDriver.hpp"
#include "AsNmpcUtils.hpp"

namespace Ipopt
{

  class AsNmpController : public AlgorithmStrategyObject
  {
    /** This is the interface for the actual controller. It handles 
     *  Data input to the controller (measurement) and returns controls */

  public: 
    
    AsNmpController(std::vector< SmartPtr<SchurDriver> >& driver_vec,
		    SmartPtr<SensitivityStepCalculator> sens_step_calc,
		    SmartPtr<Measurement> measurement,
		    Index n_nmpc_steps);

    virtual ~AsNmpController();

    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Main loop: Wait for new measurement, Get new step, maybe deal with 
     *  bounds,  see to it that everything happens in the required 
     *  timeframe. */
    NmpControllerExitStatus Run();

  private:

    std::vector< SmartPtr<SchurDriver> > driver_vec_;
    SmartPtr<SensitivityStepCalculator> sens_step_calc_;
    SmartPtr<Measurement> measurement_;
    Index n_nmpc_steps_; // I think it is useful to state this number explicitly in the constructor and here.

  };
}

#endif
