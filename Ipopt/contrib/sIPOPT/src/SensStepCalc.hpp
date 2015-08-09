// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-14

#ifndef __ASSENSSTEPCALC_HPP__
#define __ASSENSSTEPCALC_HPP__

#include "IpAlgStrategy.hpp"
#include "SensSchurDriver.hpp"


namespace Ipopt
{
  /** Forward declarations */
  class DenseVector;
  class IteratesVector;

  class SensitivityStepCalculator : public AlgorithmStrategyObject
  {
    /* This is the interface for the classes that perform the actual step. */

  public:
    SensitivityStepCalculator()
      :
      driver_(NULL),
      do_boundcheck_(false)
    {
    }

    virtual ~SensitivityStepCalculator()
    {
    }

    virtual bool InitializeImpl(const OptionsList& options,
				const std::string& prefix)
    {
      options.GetBoolValue("sens_boundcheck", do_boundcheck_, prefix);
      return true;
    }

    bool Do_Boundcheck() const
    {
      return do_boundcheck_;
    }

    void SetSchurDriver(SmartPtr<SchurDriver> driver)
    {
      DBG_ASSERT(IsValid(driver));
      driver_ = driver;
      if (IsValid(driver_->pcalc_nonconst())) {
	driver_->pcalc_nonconst()->reset_data_A();
	// when the schurdriver is set, the data in the pcalculator has to be reset to its data?
      }
    }

    SmartPtr<SchurDriver> Driver() // this should be const or protected
    {
      DBG_ASSERT(IsValid(driver_));
      return driver_;
    }

    /** This is the main algorithmic function of this class; It calculates
     *  a step using its SchurDriver, checks bounds, and returns it */
    virtual bool Step(DenseVector& delta_u, IteratesVector& sol) =0;

    /** return the sensitivity vector */
    virtual SmartPtr<IteratesVector> GetSensitivityVector() = 0;

  private:
    SmartPtr<SchurDriver> driver_;
    bool do_boundcheck_;
  };
}

#endif
