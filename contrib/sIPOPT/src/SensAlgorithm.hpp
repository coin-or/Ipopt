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

/** This is the interface for the actual controller.
 *
 *  It handles Data input to the controller (measurement) and returns controls.
 */
class SIPOPTLIB_EXPORT SensAlgorithm : public AlgorithmStrategyObject
{
public:
   SensAlgorithm(
      std::vector< SmartPtr<SchurDriver> >& driver_vec,
      SmartPtr<SensitivityStepCalculator>   sens_step_calc,
      SmartPtr<Measurement>                 measurement,
      Index                                 n_sens_steps
   );

   virtual ~SensAlgorithm();

   virtual bool InitializeImpl(
      const OptionsList& options,
      const std::string& prefix
   );

   /** Main loop: Wait for new measurement, Get new step, maybe deal with
    *  bounds,  see to it that everything happens in the required
    *  timeframe. */
   SensAlgorithmExitStatus Run();
   SensAlgorithmExitStatus ComputeSensitivityMatrix(void);

   /** accessor methods to get access to variable sizes */
   Index nl(void)
   {
      return nl_;
   }
   Index nx(void)
   {
      return nx_;
   }
   Index nzl(void)
   {
      return nzl_;
   }
   Index nzu(void)
   {
      return nzu_;
   }
   Index ns(void)
   {
      return ns_;
   }
   Index np(void)
   {
      return np_;
   }

   /** array place holders to store the vector of sensitivities */
   Number* DirectionalD_X_;
   Number* DirectionalD_L_;
   Number* DirectionalD_Z_L_;
   Number* DirectionalD_Z_U_;

   /** array place holders for the sensitivity matrix */
   Number* SensitivityM_X_;
   Number* SensitivityM_L_;
   Number* SensitivityM_Z_L_;
   Number* SensitivityM_Z_U_;

private:
   Index nl_;
   Index nx_;
   Index nzl_;
   Index nzu_;
   Index nceq_;
   Index ncineq_;
   Index ns_;
   Index np_;

   std::vector< SmartPtr<SchurDriver> > driver_vec_;
   SmartPtr<SensitivityStepCalculator> sens_step_calc_;
   SmartPtr<Measurement> measurement_;
   Index n_sens_steps_; // I think it is useful to state this number explicitly in the constructor and here.

   /** method to extract sensitivity vectors */
   void GetDirectionalDerivatives(void);

   /** method to extract sensitivity matrix */
   void GetSensitivityMatrix(
      Index col
   );

   /** private method used to uncale perturbed solution and sensitivities */
   void UnScaleIteratesVector(
      SmartPtr<IteratesVector>* V
   );
};
}

#endif
