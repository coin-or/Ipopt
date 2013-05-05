// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-06

#ifndef __SENSAPPLICATION_HPP__
#define __SENSAPPLICATION_HPP__

#include "IpReferenced.hpp"
#include "SensUtils.hpp"
#include "SensUtils.hpp"
#include "IpRegOptions.hpp"

#include "IpIpoptApplication.hpp"
#include "IpPDSystemSolver.hpp"
namespace Ipopt
{
  /** Standard exception for wrong/inconsistent suffixes for sipopt */
  DECLARE_STD_EXCEPTION(SENS_SUFFIX_ERROR);

  class SensApplication : public ReferencedObject
  {
  public:
    // constructor
    SensApplication(SmartPtr<Journalist> jnlst,
		    SmartPtr<OptionsList> options,
		    SmartPtr<RegisteredOptions> reg_options);

    ~SensApplication();

    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);

    SensAlgorithmExitStatus Run();

    void Initialize();

    void SetIpoptAlgorithmObjects(SmartPtr<IpoptApplication> app_ipopt,
				  ApplicationReturnStatus ipopt_retval);

    SmartPtr<Journalist> Jnlst()
    {
      return jnlst_;
    }


    SmartPtr<OptionsList> Options()
    {
      return options_;
    }

    /** Get the options list for setting options (const version) */
    SmartPtr<const OptionsList> Options() const
    {
      return ConstPtr(options_);
    }


  private:

    // standard constructor just so it can't be used
    //    SensApplication();

    // Pointers that are immediately passed from Ipopt and initialized by the constructor
    SmartPtr<Journalist> jnlst_;
    SmartPtr<OptionsList> options_;
    SmartPtr<IpoptData> ip_data_;
    SmartPtr<IpoptCalculatedQuantities> ip_cq_;
    SmartPtr<PDSystemSolver> pd_solver_;
    SmartPtr<IpoptNLP> ip_nlp_;
    SmartPtr<RegisteredOptions> reg_options_;
    ApplicationReturnStatus ipopt_retval_;

    /** storing options values */
    bool run_sens_;
    bool compute_red_hessian_;
    Index n_sens_steps_;
  };
}

#endif
