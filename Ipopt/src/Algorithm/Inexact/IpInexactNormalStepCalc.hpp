// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2008-08-31

#ifndef __IPINEXACTNORMALSTEPCALC_HPP__
#define __IPINEXACTNORMALSTEPCALC_HPP__

#include "IpAlgStrategy.hpp"
#include "IpInexactCq.hpp"

namespace Ipopt
{
  /** Base class for computing the normal step for the inexact step
   *  calculation algorithm.
   */
  class InexactNormalStepCalculator : public AlgorithmStrategyObject
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default onstructor */
    InexactNormalStepCalculator()
    {}

    /** Default destructor */
    virtual ~InexactNormalStepCalculator()
    {}
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix) = 0;

    /** Method for computing the normal step.  The computed step is
     *  returned as normal_x and normal_s, for the x and s variables,
     *  respectively.  These quantities are not slack-scaled.  If the
     *  step cannot be computed, this method returns false.  */
    virtual bool ComputeNormalStep(SmartPtr<Vector>& normal_x,
                                   SmartPtr<Vector>& normal_s) = 0;

  protected:
    /** Method to easily access Inexact data */
    InexactData& InexData()
    {
      InexactData& inexact_data =
        static_cast<InexactData&>(IpData().AdditionalData());
      DBG_ASSERT(dynamic_cast<InexactData*>(&IpData().AdditionalData()));
      return inexact_data;
    }

    /** Method to easily access Inexact calculated quantities */
    InexactCq& InexCq()
    {
      InexactCq& inexact_cq =
        static_cast<InexactCq&>(IpCq().AdditionalCq());
      DBG_ASSERT(dynamic_cast<InexactCq*>(&IpCq().AdditionalCq()));
      return inexact_cq;
    }

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Copy Constructor */
    InexactNormalStepCalculator(const InexactNormalStepCalculator&);

    /** Overloaded Equals Operator */
    void operator=(const InexactNormalStepCalculator&);
    //@}
  };

} // namespace Ipopt

#endif
