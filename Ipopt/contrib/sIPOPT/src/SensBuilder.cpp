// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-10

#include "AsSchurBuilder.hpp"
#include "AsPCalculator.hpp"
#include "AsStdPCalculator.hpp"
#include "AsIndexPCalculator.hpp"
#include "AsSchurData.hpp"
#include "AsStdSchurData.hpp"
#include "AsIndexSchurData.hpp"
#include "AsDenseGenSchurDriver.hpp"
#include "AsIFTSchurDriver.hpp"
#include "AsMeasurement.hpp"
#include "AsMetadataMeasurement.hpp"
#include "AsStdStepCalc.hpp"

#include <string>
#include <sstream>

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 1;
#endif

  SchurBuilder::SchurBuilder()
  {
    DBG_START_METH("SchurBuilder::SchurBuilder", dbg_verbosity);
  }

  SchurBuilder::~SchurBuilder()
  {
    DBG_START_METH("SchurBuilder::~SchurBuilder", dbg_verbosity);
  }

  SmartPtr<AsNmpController> SchurBuilder::BuildNmpc(const Journalist& jnlst,
						    const OptionsList& options,
						    const std::string& prefix,
						    IpoptNLP& ip_nlp,
						    IpoptData& ip_data,
						    IpoptCalculatedQuantities& ip_cq,
						    PDSystemSolver& pd_solver)
  {
    DBG_START_METH("SchurBuilder::BuildNmpc", dbg_verbosity);

    // Very first thing is setting trial = curr.
    SmartPtr<IteratesVector> trialcopyvector = ip_data.curr()->MakeNewIteratesVectorCopy();
    ip_data.set_trial(trialcopyvector);

    // Check options which Backsolver to use here
    SmartPtr<AsBacksolver> backsolver = new SimpleBacksolver(&pd_solver);

    // Check option which SchurData and PCalculator implementation to use
    std::string nmpc_calc_style;
    options.GetStringValue("nmpc_calc_style", nmpc_calc_style, "");

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
    if (nmpc_calc_style=="index") {
      E_0 = new IndexSchurData();
    }
    else if (nmpc_calc_style=="std") {
      E_0 = new StdSchurData();
    }
    
    std::string select_step;
    options.GetStringValue("select_step", select_step, "");
    DBG_PRINT((dbg_verbosity, "Selected step is %s\n", select_step.c_str()));

    std::vector<Index> initial_c = measurement->GetInitialEqConstraints(); // type: List
    if (select_step=="advanced") {
      std::vector<Index> z_k_index = measurement->GetNmpcState(1); // type: Index 
    
      E_0->SetData_Index(z_k_index.size(), &z_k_index[0]);
    
      std::vector<Index> delta_u_sort_empty;
      Index new_du_size_empty=0;
    
      E_0->AddData_List(initial_c, delta_u_sort_empty,new_du_size_empty,1);
    }
    else if (select_step=="sensitivity" || select_step=="iftsensitivity") {
      E_0->SetData_List(initial_c);
    }
    else if (select_step=="ift") {
      /* For the IFT-step, the variables have to be ordered the other way around:
	 The delta_u for the lambda_0 goes into the line of z_0 
	 and vice versa. */
      std::vector<Index> z_k_index = measurement->GetNmpcState(1); // type: Index 
      Index i_c_size = initial_c.size();
      Index count = i_c_size;
      std::vector<Index>::iterator it=z_k_index.begin();
      while (it!=z_k_index.end() && count<2*i_c_size) {
	if (*it>0) {
	  initial_c.push_back(*it+i_c_size);
	}
	++it;
      }
      E_0->SetData_List(initial_c);
    }

    E_0->Print(jnlst,J_VECTOR,J_USER1,"E_0");

    SmartPtr<PCalculator> pcalc;
    if (select_step=="advanced" || select_step=="sensitivity") { // don't create pcalculator for IFT step
      // Check options which PCalculator to use here
      if (nmpc_calc_style=="index") {
	pcalc = new IndexPCalculator(backsolver, E_0);
      }
      else if (nmpc_calc_style=="std") {
	pcalc = new StdPCalculator(backsolver, E_0);
      }

      bool retval = pcalc->Initialize(jnlst,
				      ip_nlp,
				      ip_data,
				      ip_cq,
				      options,
				      prefix);
      DBG_ASSERT(retval);
      retval = pcalc->ComputeP();
      //pcalc->Print(jnlst,J_VECTOR,J_USER1,"PCalc");
      if (!retval) {
	// Throw exception that P calculation failed
      }
    }
    else { // ift steps get an empty pcalculator if boundcheck is on (for fix-relax)
      bool bound_check;
      options.GetBoolValue("sens_boundcheck", bound_check, prefix);
      if (bound_check) {
	if (nmpc_calc_style=="index") {
	  pcalc = new IndexPCalculator(backsolver, new IndexSchurData());
	  bool retval = pcalc->Initialize(jnlst,
					  ip_nlp,
					  ip_data,
					  ip_cq,
					  options,
					  prefix);
	  DBG_ASSERT(retval);
      	}
      }
    }
    

    // Find out how many steps there are and create as many SchurSolveDrivers
    int n_sens_steps;
    options.GetIntegerValue("n_sens_steps",n_sens_steps,prefix);
    //DBG_ASSERT(n_sens_steps<2); // for testing the new formula, can't do more!

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
      if (select_step=="advanced") {
	assert(false);
	//driver_vec[i] = new DenseGenSchurDriver(backsolver, pcalc, E_0);
      }
      else if (select_step=="ift" || "iftsensitivity") {
	driver_vec[i] = new IFTSchurDriver(backsolver, pcalc,E_0);
      } else if (select_step=="sensitivity") {
	// Create SchurDriver from pcalc and suffix indices
	SmartPtr<SchurData> E_i;
	if (nmpc_calc_style=="index") { 
	  E_i = new IndexSchurData();
	}
	else if (nmpc_calc_style=="std") {
	  E_i = new StdSchurData();
	}
	sens_state_list = measurement->GetNmpcState(i+1);
	DBG_PRINT((dbg_verbosity, "sens_state_list.size()=%d", sens_state_list.size()));
	// E_i->SetData_List(sens_state_list); // this is obsolete since Measurement class changed behaviour and now outputs indices!
	E_i->SetData_Index(sens_state_list.size(),&sens_state_list[0]);
	E_i_name = "E_";
	append_Index(E_i_name, i+1);
	E_i->Print(jnlst,J_VECTOR,J_USER1,E_i_name.c_str());
	assert(false);
	//driver_vec[i] = new DenseGenSchurDriver(backsolver, pcalc, E_i);
      }
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

    SmartPtr<AsNmpController> controller = new AsNmpController(driver_vec,
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

  SmartPtr<ReducedHessianCalculator> SchurBuilder::BuildRedHessCalc(const Journalist& jnlst,
								    const OptionsList& options,
								    const std::string& prefix,
								    IpoptNLP& ip_nlp,
								    IpoptData& ip_data,
								    IpoptCalculatedQuantities& ip_cq,
								    //NmpcTNLPAdapter& nmpc_tnlp_adapter,
								    PDSystemSolver& pd_solver)
  {
    DBG_START_METH("SchurBuilder::BuildRedHessCalc", dbg_verbosity);
    
    // Check options which Backsolver to use here
    SmartPtr<AsBacksolver> backsolver = new SimpleBacksolver(&pd_solver);

    // Check option which SchurData and PCalculator implementation to use
    std::string nmpc_calc_style;
    options.GetStringValue("nmpc_calc_style", nmpc_calc_style, "");

    // Create suffix handler
    SmartPtr<SuffixHandler> suffix_handler = new MetadataMeasurement();
    dynamic_cast<MetadataMeasurement*>(GetRawPtr(suffix_handler))->Initialize(jnlst,
									      ip_nlp,
									      ip_data,
									      ip_cq,
									      options,
									      prefix);
    SmartPtr<SchurData> E_0;
    if (nmpc_calc_style=="index") {
      E_0 = new IndexSchurData();
    }
    else if (nmpc_calc_style=="std") {
      E_0 = new StdSchurData();
    }
    
    std::vector<Index> hessian_suff = suffix_handler->GetIntegerSuffix("red_hessian");

    Index setdata_error = E_0->SetData_Index(hessian_suff.size(), &hessian_suff[0], 1.0);
    if ( setdata_error ){
      jnlst.Printf(J_ERROR, J_MAIN, "\nEXIT: An Error Occured while processing "
		    "the Indices for the reduced hessian computation: Something "
		    "is wrong with index %d\n",setdata_error);
      THROW_EXCEPTION(ASNMPC_BUILDER_ERROR,
                      "Reduced Hessian Index Error");
    }
    
    SmartPtr<PCalculator> pcalc;
    if (nmpc_calc_style=="index") {
      pcalc = new IndexPCalculator(backsolver, E_0);
    }
    else if (nmpc_calc_style=="std") {
      pcalc = new StdPCalculator(backsolver, E_0);
    }

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

