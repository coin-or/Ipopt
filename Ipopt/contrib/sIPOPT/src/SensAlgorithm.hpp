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


    /** accessor methods to get access to variable sizes */
    Index nl(void) { return nl_ ; }
    Index nx(void) { return nx_ ; }
    Index nzl(void) {return nzl_ ; }
    Index nzu(void) {return nzu_ ; }
    
    /** array place holders to store the vector of sensitivities */
    Number *Sensitivity_X_ ;
    Number *Sensitivity_L_ ;
    Number *Sensitivity_Z_U_ ;
    Number *Sensitivity_Z_L_ ;

  private:
    
    Index nl_ ;
    Index nx_ ;  
    Index nzl_ ; 
    Index nzu_ ; 
    Index nceq_ ;
    Index ncineq_ ;

    std::vector< SmartPtr<SchurDriver> > driver_vec_;
    SmartPtr<SensitivityStepCalculator> sens_step_calc_;
    SmartPtr<Measurement> measurement_;
    Index n_sens_steps_; // I think it is useful to state this number explicitly in the constructor and here.

    /** method to extract sensitivity vectors */
    void GetSensitivities(void) ;

    /** private method used to uncale perturbed solution and sensitivities */
    void UnScaleIteratesVector(SmartPtr<IteratesVector> *V) ;
  };
}

#endif
