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
#include "SensAlgorithm.hpp"
#include "IpRegOptions.hpp"

#include "IpIpoptApplication.hpp"
#include "IpPDSystemSolver.hpp"

#include "IpSmartPtr.hpp"

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

    /** Copy over value of ds/dp */
    void GetSensitivities(Number *SX, Number *SL, Number *SZL, Number *SZU) {
      if (GetRawPtr(controller) != NULL && NULL != Sensitivity_X && NULL != Sensitivity_Z_L && NULL != Sensitivity_Z_U && NULL != Sensitivity_L) { 
	
	for (int i=0 ; i < controller->nx(); ++i) SX[i] = Sensitivity_X[i] ;
	for (int i=0; i < controller->nzu(); ++i) SZU[i] = Sensitivity_Z_U[i] ;
	for (int i=0; i < controller->nzl(); ++i) SZL[i] = Sensitivity_Z_L[i] ;
	for (int i=0 ; i < controller->nl(); ++i) SL[i]  = Sensitivity_L[i] ;
      
      }
    }

    /** accessor methods to get sizing info */
    Index nx()  {return (GetRawPtr(controller)!=NULL) ? controller->nx() : -1 ;} 
    Index nl()  {return (GetRawPtr(controller)!=NULL) ? controller->nl() : -1 ;} 
    Index nzu() {return (GetRawPtr(controller)!=NULL) ? controller->nzu(): -1 ;} 
    Index nzl() {return (GetRawPtr(controller)!=NULL) ? controller->nzl(): -1 ;} 

    /* place holders to keep the values of ds/dp for each type of variable */
    Number *Sensitivity_X ;
    Number *Sensitivity_L ;
    Number *Sensitivity_Z_U ;
    Number *Sensitivity_Z_L ;

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

    SmartPtr<SensAlgorithm> controller ;

    /** storing options values */
    bool run_sens_;
    bool compute_red_hessian_;
    Index n_sens_steps_;
  };
}

#endif
