// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2008-08-31

#ifndef __IPINEXACTDOGLEGNORMAL_HPP__
#define __IPINEXACTDOGLEGNORMAL_HPP__

#include "IpInexactNormalStepCalc.hpp"
#include "IpInexactNewtonNormal.hpp"

namespace Ipopt
{
  /** Compute the normal step using a dogleg approach.
   */
  class InexactDoglegNormalStep : public InexactNormalStepCalculator
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default onstructor */
    InexactDoglegNormalStep(SmartPtr<InexactNewtonNormalStep> newton_step);

    /** Default destructor */
    virtual ~InexactDoglegNormalStep();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Method for computing the normal step.  The computed step is
     *  returned as normal_x and normal_s, for the x and s variables,
     *  respectively.  These quantities are not slack-scaled.  If the
     *  step cannot be computed, this method returns false.  */
    virtual bool ComputeNormalStep(SmartPtr<Vector>& normal_x,
                                   SmartPtr<Vector>& normal_s);

    /** Methods for IpoptType */
    //@{
    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);
    //@}

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
    InexactDoglegNormalStep();

    /** Copy Constructor */
    InexactDoglegNormalStep(const InexactDoglegNormalStep&);

    /** Overloaded Equals Operator */
    void operator=(const InexactDoglegNormalStep&);
    //@}

    /** Point to object for computing the "Newton" step in the dogleg
     *  method */
    SmartPtr<InexactNewtonNormalStep> newton_step_;

    /** @name Algorithmic options */
    //@{
    Number omega_max_;
    //@}

    /** Current value of the trust region factor */
    Number curr_omega_;

    /** Flag indicating if trust region was active in last iteration */
    bool last_tr_inactive_;
  };

} // namespace Ipopt

#endif
