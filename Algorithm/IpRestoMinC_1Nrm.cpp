// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpRestoMinC_1Nrm.hpp"
#include "IpCompoundVector.hpp"
#include "IpRestoIpoptNLP.hpp"

namespace Ipopt
{

  MinC_1NrmRestorationPhase::MinC_1NrmRestorationPhase
  (IpoptAlgorithm& resto_alg)
      :
      resto_alg_(&resto_alg),
      resto_options_(NULL)
  {
    DBG_ASSERT(IsValid(resto_alg_));
  }

  MinC_1NrmRestorationPhase::~MinC_1NrmRestorationPhase()
  {}

  bool MinC_1NrmRestorationPhase::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    // keep a copy of these options to use when setting up the
    // restoration phase
    resto_options_ = new OptionsList(options);

    return true;
  }

  bool
  MinC_1NrmRestorationPhase::PerformRestoration()
  {
    Jnlst().Printf(J_DETAILED, J_MAIN, "Starting Restoration Phase\n");

    // Create the restoration phase NLP etc objects
    SmartPtr<IpoptData> resto_ip_data = new IpoptData();
    SmartPtr<IpoptNLP> resto_ip_nlp
    = new RestoIpoptNLP(IpNLP(), IpData(), IpCq(), *resto_ip_data);
    SmartPtr<IpoptCalculatedQuantities> resto_ip_cq
    = new IpoptCalculatedQuantities(resto_ip_nlp, resto_ip_data);

    // Initialize the restoration phase algorithm
    resto_alg_->Initialize(Jnlst(), *resto_ip_nlp, *resto_ip_data,
                           *resto_ip_cq, *resto_options_, "resto.");

    // Set iteration counter and info field for the restoration phase
    resto_ip_data->Set_iter_count(IpData().iter_count()+1);
    resto_ip_data->Set_info_regu_x(IpData().info_regu_x());
    resto_ip_data->Set_info_alpha_primal(IpData().info_alpha_primal());
    resto_ip_data->Set_info_alpha_primal_char(IpData().info_alpha_primal_char());
    resto_ip_data->Set_info_alpha_dual(IpData().info_alpha_dual());
    resto_ip_data->Set_info_ls_count(IpData().info_ls_count());

    IpoptAlgorithm::SolverReturn resto_status	= resto_alg_->Optimize();

    int retval=-1;

    if (resto_status == IpoptAlgorithm::SUCCESS) {
      if (Jnlst().ProduceOutput(J_DETAILED, J_SOLUTION)) {
        Jnlst().Printf(J_DETAILED, J_SOLUTION,
                       "\nRESTORATION PHASE RESULTS\n");
        Jnlst().Printf(J_DETAILED, J_SOLUTION,
                       "\n\nOptimal solution found! \n");
        Jnlst().Printf(J_DETAILED, J_SOLUTION,
                       "Optimal Objective Value = %.16E\n", resto_ip_cq->curr_f());
        Jnlst().Printf(J_DETAILED, J_SOLUTION,
                       "Number of Iterations = %d\n", resto_ip_data->iter_count());
      }
      if (Jnlst().ProduceOutput(J_VECTOR, J_SOLUTION)) {
        Jnlst().PrintVector(J_VECTOR, J_SOLUTION,
                            "x", *resto_ip_data->curr_x());
        Jnlst().PrintVector(J_VECTOR, J_SOLUTION,
                            "y_c", *resto_ip_data->curr_y_c());
        Jnlst().PrintVector(J_VECTOR, J_SOLUTION,
                            "y_d", *resto_ip_data->curr_y_d());
        Jnlst().PrintVector(J_VECTOR, J_SOLUTION,
                            "z_L", *resto_ip_data->curr_z_L());
        Jnlst().PrintVector(J_VECTOR, J_SOLUTION,
                            "z_U", *resto_ip_data->curr_z_U());
        Jnlst().PrintVector(J_VECTOR, J_SOLUTION,
                            "v_L", *resto_ip_data->curr_v_L());
        Jnlst().PrintVector(J_VECTOR, J_SOLUTION,
                            "v_U", *resto_ip_data->curr_v_U());
      }

      retval = 0;
    }
    else {
      Jnlst().Printf(J_ERROR, J_MAIN, "Sorry, things failed ?!?!\n");
      retval = 1;
    }

    if (retval == 0) {
      // Copy the results back...!!!
      // for now, I only copy x and s
      SmartPtr<const CompoundVector> cx =
        dynamic_cast<const CompoundVector*>(GetRawPtr(resto_ip_data->curr_x()));
      DBG_ASSERT(IsValid(cx));
      SmartPtr<const Vector> x_only = cx->GetComp(0);
      SmartPtr<const Vector> s_only = resto_ip_data->curr_s();
      IpData().SetTrialPrimalVariablesFromPtr(x_only, s_only);

      IpData().Set_iter_count(resto_ip_data->iter_count()-1);
      // Skip the next line, because it would just replicate the first
      // on during the restoration phase.
      IpData().Set_info_skip_output(true);
    }

    return (retval == 0);
  }

} // namespace Ipopt
