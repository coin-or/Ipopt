// Copyright (C) 2008, 2011 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2008-09-05

#ifndef __IPINEXACTNEWTONNORMAL_HPP__
#define __IPINEXACTNEWTONNORMAL_HPP__

#include "IpAlgStrategy.hpp"
#include "IpAugSystemSolver.hpp"
#include "IpInexactCq.hpp"

namespace Ipopt
{
  /** Compute the "Newton" normal step from the (slack-scaled)
   *  augmented system.
   */
  class InexactNewtonNormalStep : public AlgorithmStrategyObject
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default onstructor */
    InexactNewtonNormalStep(SmartPtr<AugSystemSolver> aug_solver);

    /** Default destructor */
    virtual ~InexactNewtonNormalStep();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Method for computing the normal step.  The computed step is
     *  returned as normal_x and normal_s, for the x and s variables,
     *  respectively.  These quantities are not in the original space,
     *  but in the space scaled by the slacks.  If the step cannot be
     *  computed, this method returns false.  */
    virtual bool ComputeNewtonNormalStep(Vector& newton_x, Vector& newton_s);

    /** Methods for IpoptType */
    //@{
    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);
    //@}

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
    /** Default onstructor */
    InexactNewtonNormalStep();

    /** Copy Constructor */
    InexactNewtonNormalStep(const InexactNewtonNormalStep&);

    /** Overloaded Equals Operator */
    void operator=(const InexactNewtonNormalStep&);
    //@}

    /** Object to be used to solve the augmented system */
    SmartPtr<AugSystemSolver> aug_solver_;
  };

} // namespace Ipopt

#endif
