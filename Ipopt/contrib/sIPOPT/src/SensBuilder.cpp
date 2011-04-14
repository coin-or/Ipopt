// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-10

#include "SensBuilder.hpp"
#include "SensPCalculator.hpp"
#include "SensIndexPCalculator.hpp"
#include "SensSchurData.hpp"
#include "SensIndexSchurData.hpp"
#include "SensDenseGenSchurDriver.hpp"
#include "SensMeasurement.hpp"
#include "SensMetadataMeasurement.hpp"
#include "SensStdStepCalc.hpp"

#include <string>
#include <sstream>

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 1;
#endif

  SensBuilder::SensBuilder()
  {
    DBG_START_METH("SensBuilder::SensBuilder", dbg_verbosity);
  }

  SensBuilder::~SensBuilder()
  {
    DBG_START_METH("SensBuilder::~SensBuilder", dbg_verbosity);
  }

  SmartPtr<SensAlgorithm> SensBuilder::BuildSensAlg(const Journalist& jnlst,
						    const OptionsList& options,
						    const std::string& prefix,
						    IpoptNLP& ip_nlp,
						    IpoptData& ip_data,
						    IpoptCalculatedQuantities& ip_cq,
						    PDSystemSolver& pd_solver)
  {
    DBG_START_METH("SensBuilder::BuildSensAlg", dbg_verbosity);

    // Very first thing is setting trial = curr.
    SmartPtr<IteratesVector> trialcopyvector = ip_data.curr()->MakeNewIteratesVectorCopy();
    ip_data.set_trial(trialcopyvector);

    // Check options which Backsolver to use here
    SmartPtr<SensBacksolver> backsolver = new SimpleBacksolver(&pd_solver);

    // Create measurement unit
    SmartPtr<Measurement> measurement = new MetadataMeasurement();
    (dynamic_cast<MetadataMeasurement*>(GetRawPtr(measurement)))->Initialize(jnlst,
									     ip_nlp,
									     ip_data,
									     ip_cq,
									     options,
									     prefix);

    // Check ParameterData, send it to Pcalculator
    SmartPtr<SchurData> E_0;
    E_0 = new IndexSchurData();

    std::vector<Index> initial_c = measurement->GetInitialEqConstraints(); // type: List
    E_0->SetData_List(initial_c);
    E_0->Print(jnlst,J_VECTOR,J_USER1,"E_0");

    SmartPtr<PCalculator> pcalc;
    bool bound_check;
    options.GetBoolValue("sens_boundcheck", bound_check, prefix);
    if (bound_check) {
      pcalc = new IndexPCalculator(backsolver, new IndexSchurData());
      bool retval = pcalc->Initialize(jnlst,
				      ip_nlp,
				      ip_data,
				      ip_cq,
				      options,
				      prefix);
      DBG_ASSERT(retval);
    }

    // Find out how many steps there are and create as many SchurSolveDrivers
    int n_sens_steps;
    options.GetIntegerValue("n_sens_steps",n_sens_steps,prefix);

    // Create std::vector container in which we are going to keep the SchurDrivers
    std::vector< SmartPtr<SchurDriver> > driver_vec(n_sens_steps);

    /** Here there should be the point to pass on the driver_vec and fork off the
     *  Schurcomputations to a different function/process if needed. */
    std::vector<Index> sens_state_list;
    Index schur_retval;
    std::string E_i_name;

    /** THIS FOR-LOOP should be done better with a better
     *  Measurement class. This should get it's own branch! */
    for (Index i=0; i<n_sens_steps; ++i) {
      driver_vec[i] = new DenseGenSchurDriver(backsolver, pcalc,E_0);
      driver_vec[i]->Initialize(jnlst,
				ip_nlp,
				ip_data,
				ip_cq,
				options,
				prefix);
      schur_retval = driver_vec[i]->SchurBuild();
      DBG_ASSERT(schur_retval);
      schur_retval = driver_vec[i]->SchurFactorize();
      DBG_ASSERT(schur_retval);
    }

    SmartPtr<SensitivityStepCalculator> sens_stepper = new StdStepCalculator(E_0, backsolver);

    sens_stepper->Initialize(jnlst,
			     ip_nlp,
			     ip_data,
			     ip_cq,
			     options,
			     prefix);

    SmartPtr<SensAlgorithm> controller = new SensAlgorithm(driver_vec,
							   sens_stepper,
							   measurement,
							   n_sens_steps);

    controller->Initialize(jnlst,
			   ip_nlp,
			   ip_data,
			   ip_cq,
			   options,
			   prefix);
    return controller;
  }

  SmartPtr<ReducedHessianCalculator> SensBuilder::BuildRedHessCalc(const Journalist& jnlst,
								   const OptionsList& options,
								   const std::string& prefix,
								   IpoptNLP& ip_nlp,
								   IpoptData& ip_data,
								   IpoptCalculatedQuantities& ip_cq,
								   PDSystemSolver& pd_solver)
  {
    DBG_START_METH("SensBuilder::BuildRedHessCalc", dbg_verbosity);

    // Check options which Backsolver to use here
    SmartPtr<SensBacksolver> backsolver = new SimpleBacksolver(&pd_solver);

    // Create suffix handler
    SmartPtr<SuffixHandler> suffix_handler = new MetadataMeasurement();
    dynamic_cast<MetadataMeasurement*>(GetRawPtr(suffix_handler))->Initialize(jnlst,
									      ip_nlp,
									      ip_data,
									      ip_cq,
									      options,
									      prefix);
    SmartPtr<SchurData> E_0;
    E_0 = new IndexSchurData();

    std::vector<Index> hessian_suff = suffix_handler->GetIntegerSuffix("red_hessian");

    Index setdata_error = E_0->SetData_Index(hessian_suff.size(), &hessian_suff[0], 1.0);
    if ( setdata_error ){
      jnlst.Printf(J_ERROR, J_MAIN, "\nEXIT: An Error Occured while processing "
		   "the Indices for the reduced hessian computation: Something "
		   "is wrong with index %d\n",setdata_error);
      THROW_EXCEPTION(SENS_BUILDER_ERROR,
                      "Reduced Hessian Index Error");
    }

    SmartPtr<PCalculator> pcalc;
    pcalc = new IndexPCalculator(backsolver, E_0);

    bool retval = pcalc->Initialize(jnlst,
				    ip_nlp,
				    ip_data,
				    ip_cq,
				    options,
				    prefix);
    DBG_ASSERT(retval);

    retval = pcalc->ComputeP();

    SmartPtr<ReducedHessianCalculator> red_hess_calc = new ReducedHessianCalculator(E_0, pcalc);

    retval = red_hess_calc->Initialize(jnlst,
				       ip_nlp,
				       ip_data,
				       ip_cq,
				       options,
				       prefix);

    return red_hess_calc;
  }
}
