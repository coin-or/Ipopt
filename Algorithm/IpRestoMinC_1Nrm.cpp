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
  static const Index dbg_verbosity = 0;

  MinC_1NrmRestorationPhase::MinC_1NrmRestorationPhase
  (IpoptAlgorithm& resto_alg,
   const SmartPtr<EqMultiplierCalculator>& eq_mult_calculator)
      :
      resto_alg_(&resto_alg),
      eq_mult_calculator_(eq_mult_calculator),
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

    Number value;
    if (options.GetNumericValue("laminitmax", value, prefix)) {
      ASSERT_EXCEPTION(value >= 0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"laminitmax\": Value must be non-negative.");
      laminitmax_ = value;
    }
    else {
      laminitmax_ = 1e3;
    }

    bool retvalue = true;
    if (IsValid(eq_mult_calculator_)) {
      retvalue = eq_mult_calculator_->Initialize(Jnlst(), IpNLP(), IpData(),
                 IpCq(), options, prefix);
    }
    return retvalue;
  }

  bool
  MinC_1NrmRestorationPhase::PerformRestoration()
  {
    DBG_START_METH("MinC_1NrmRestorationPhase::PerformRestoration",
                   dbg_verbosity);
    Jnlst().Printf(J_DETAILED, J_MAIN, "Starting Restoration Phase\n");

    DBG_ASSERT(IpCq().curr_constraint_violation()>0.);

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

      // Recompute the equality constraint multipliers as least square estimate
      if (IsValid(eq_mult_calculator_) && laminitmax_>0.) {
        // First move all the trial data into the current fields, since
        // those values are needed to compute the initial values for
        // the multipliers
        SmartPtr<Vector> y_c = IpData().curr_y_c()->MakeNew();
        SmartPtr<Vector> y_d = IpData().curr_y_d()->MakeNew();
        IpData().CopyTrialToCurrent();
        y_c = IpData().curr_y_c()->MakeNew();
        y_d = IpData().curr_y_d()->MakeNew();
        bool retval = eq_mult_calculator_->CalculateMultipliers(*y_c, *y_d);
        if (!retval) {
          y_c->Set(0.0);
          y_d->Set(0.0);
        }
        else {
          Jnlst().Printf(J_DETAILED, J_INITIALIZATION,
                         "Least square estimates max(y_c) = %e, max(y_d) = %e\n",
                         y_c->Amax(), y_d->Amax());
          Number laminitnrm = Max(y_c->Amax(), y_d->Amax());
          if (!retval || laminitnrm > laminitmax_) {
            y_c->Set(0.0);
            y_d->Set(0.0);
          }
        }
        IpData().SetTrialEqMultipliers(*y_c, *y_d);
      }
      else {
        SmartPtr<Vector> y_c = IpData().curr_y_c()->MakeNew();
        SmartPtr<Vector> y_d = IpData().curr_y_d()->MakeNew();
        y_c->Set(0.0);
        y_d->Set(0.0);
        IpData().SetTrialEqMultipliers(*y_c, *y_d);
      }

      DBG_PRINT_VECTOR(2, "y_c", *IpData().curr_y_c());
      DBG_PRINT_VECTOR(2, "y_d", *IpData().curr_y_d());

      IpData().Set_iter_count(resto_ip_data->iter_count()-1);
      // Skip the next line, because it would just replicate the first
      // on during the restoration phase.
      IpData().Set_info_skip_output(true);
    }

    return (retval == 0);
  }

} // namespace Ipopt
