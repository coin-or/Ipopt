// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-10

#ifndef __ASSCHURBUILDER_HPP__
#define __ASSCHURBUILDER_HPP__


#include "IpReferenced.hpp"
#include "SensAlgorithm.hpp"
#include "IpPDSystemSolver.hpp"
#include "SensUtils.hpp"
#include "SensReducedHessianCalculator.hpp"


namespace Ipopt
{
  DECLARE_STD_EXCEPTION(SENS_BUILDER_ERROR);

  class SensBuilder : public ReferencedObject
  {
    /** This class sets up everything necessary and
     *  builds the P matrix which is an intermediate step
     *  in calculating the schur matrix. */
  public:
    SensBuilder();

    ~SensBuilder();

    SmartPtr<SensAlgorithm> BuildSensAlg(const Journalist& jnlst,
					 const OptionsList& options,
					 const std::string& prefix,
					 IpoptNLP& ip_nlp,
					 IpoptData& ip_data,
					 IpoptCalculatedQuantities& ip_cq,
					 PDSystemSolver& pd_solver);

    SmartPtr<ReducedHessianCalculator> BuildRedHessCalc(const Journalist& jnlst,
							const OptionsList& options,
							const std::string& prefix,
							IpoptNLP& ip_nlp,
							IpoptData& ip_data,
							IpoptCalculatedQuantities& ip_cq,
							PDSystemSolver& pd_solver);

  private:

  };

}

#endif
