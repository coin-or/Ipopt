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
  DBG_SET_VERBOSITY(0);

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

    if (options.GetNumericValue("boundmultinitmax", value, prefix)) {
      ASSERT_EXCEPTION(value >= 0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"boundmultinitmax\": Value must be non-negative.");
      boundmultinitmax_ = value;
    }
    else {
      boundmultinitmax_ = 1e3;
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

    // Call the optimization algorithm to solve the restoration phase
    // problem
    IpoptAlgorithm::SolverReturn resto_status	= resto_alg_->Optimize();

    int retval=-1;

    if (resto_status == IpoptAlgorithm::SUCCESS) {
      if (Jnlst().ProduceOutput(J_DETAILED, J_LINE_SEARCH)) {
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                       "\nRESTORATION PHASE RESULTS\n");
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                       "\n\nOptimal solution found! \n");
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                       "Optimal Objective Value = %.16E\n", resto_ip_cq->curr_f());
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                       "Number of Iterations = %d\n", resto_ip_data->iter_count());
      }
      if (Jnlst().ProduceOutput(J_VECTOR, J_LINE_SEARCH)) {
        Jnlst().PrintVector(J_VECTOR, J_LINE_SEARCH,
                            "x", *resto_ip_data->curr_x());
        Jnlst().PrintVector(J_VECTOR, J_LINE_SEARCH,
                            "y_c", *resto_ip_data->curr_y_c());
        Jnlst().PrintVector(J_VECTOR, J_LINE_SEARCH,
                            "y_d", *resto_ip_data->curr_y_d());
        Jnlst().PrintVector(J_VECTOR, J_LINE_SEARCH,
                            "z_L", *resto_ip_data->curr_z_L());
        Jnlst().PrintVector(J_VECTOR, J_LINE_SEARCH,
                            "z_U", *resto_ip_data->curr_z_U());
        Jnlst().PrintVector(J_VECTOR, J_LINE_SEARCH,
                            "v_L", *resto_ip_data->curr_v_L());
        Jnlst().PrintVector(J_VECTOR, J_LINE_SEARCH,
                            "v_U", *resto_ip_data->curr_v_U());
      }

      retval = 0;
    }
    else if (resto_status == IpoptAlgorithm::MAXITER_EXCEEDED) {
      //ToDo
      THROW_EXCEPTION(IpoptException, "Maximal number of iterations exceeded in restoration phase.");
      retval = 1;
    }
    else {
      Jnlst().Printf(J_ERROR, J_MAIN, "Sorry, things failed ?!?!\n");
      retval = 1;
    }

    if (retval == 0) {
      // Copy the results into the trial fields;. They will be
      // accepted later in the full algorith
      SmartPtr<const CompoundVector> cx =
        dynamic_cast<const CompoundVector*>(GetRawPtr(resto_ip_data->curr_x()));
      DBG_ASSERT(IsValid(cx));
      SmartPtr<const Vector> x_only = cx->GetComp(0);
      SmartPtr<const Vector> s_only = resto_ip_data->curr_s();
      IpData().SetTrialPrimalVariablesFromPtr(x_only, s_only);

      // Update the bound multiplers, pretending that the entire
      // progress in x and s in the restoration phase has been one
      // [rimal-dual Newton step (and therefore the result of solving
      // an augmented system)
      SmartPtr<Vector> delta_z_L = IpData().curr_z_L()->MakeNew();
      SmartPtr<Vector> delta_z_U = IpData().curr_z_U()->MakeNew();
      SmartPtr<Vector> delta_v_L = IpData().curr_v_L()->MakeNew();
      SmartPtr<Vector> delta_v_U = IpData().curr_v_U()->MakeNew();
      ComputeBoundMultiplierStep(*delta_z_L, *IpData().curr_z_L(),
				 *IpCq().curr_slack_x_L(),
				 *IpCq().trial_slack_x_L());
      ComputeBoundMultiplierStep(*delta_z_U, *IpData().curr_z_U(),
				 *IpCq().curr_slack_x_U(),
				 *IpCq().trial_slack_x_U());
      ComputeBoundMultiplierStep(*delta_v_L, *IpData().curr_v_L(),
				 *IpCq().curr_slack_s_L(),
				 *IpCq().trial_slack_s_L());
      ComputeBoundMultiplierStep(*delta_v_U, *IpData().curr_v_U(),
				 *IpCq().curr_slack_s_U(),
				 *IpCq().trial_slack_s_U());

      DBG_PRINT_VECTOR(1, "delta_z_L", *delta_z_L);
      DBG_PRINT_VECTOR(1, "delta_z_U", *delta_z_U);
      DBG_PRINT_VECTOR(1, "delta_v_L", *delta_v_L);
      DBG_PRINT_VECTOR(1, "delta_v_U", *delta_v_U);

      Number alpha_dual = IpCq().dual_frac_to_the_bound(IpData().curr_tau(),
							*delta_z_L,
							*delta_z_U,
							*delta_v_L,
							*delta_v_U);
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
		     "Step size for bound multipliers: %8.2e\n", alpha_dual);
     
      IpData().SetTrialBoundMultipliersFromStep(alpha_dual,
						*delta_z_L,
						*delta_z_U,
						*delta_v_L,
						*delta_v_U);
      
#ifdef olddd
      // DELETEME
      // ToDo: For the bound multipliers, for now we just keep the
      // current ones
      IpData().SetTrialBoundMultipliersFromPtr(IpData().curr_z_L(),
          IpData().curr_z_U(),
          IpData().curr_v_L(),
          IpData().curr_v_U());
#endif

      // ToDo: Check what to do here:
      Number boundmultmax = Max(IpData().trial_z_L()->Amax(),
				IpData().trial_z_U()->Amax(),
				IpData().trial_v_L()->Amax(),
				IpData().trial_v_U()->Amax());
      if (boundmultmax > boundmultinitmax_) {
	Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
		       "Bound multipliers after restoration phase too large (max=%8.2e). Set all to 1.\n",
		       boundmultmax);
	SmartPtr<Vector> new_z_L = IpData().curr_z_L()->MakeNew();
	SmartPtr<Vector> new_z_U = IpData().curr_z_U()->MakeNew();
	SmartPtr<Vector> new_v_L = IpData().curr_v_L()->MakeNew();
	SmartPtr<Vector> new_v_U = IpData().curr_v_U()->MakeNew();
	new_z_L->Set(1.0);
	new_z_U->Set(1.0);
	new_v_L->Set(1.0);
	new_v_U->Set(1.0);      
	IpData().SetTrialBoundMultipliersFromPtr(GetRawPtr(new_z_L), GetRawPtr(new_z_U),
						 GetRawPtr(new_v_L), GetRawPtr(new_v_U));
	
      }
      // Recompute the equality constraint multipliers as least square estimate
      if (IsValid(eq_mult_calculator_) && laminitmax_>0.) {
        // First move all the trial data into the current fields, since
        // those values are needed to compute the initial values for
        // the multipliers
        IpData().CopyTrialToCurrent();
        SmartPtr<Vector> y_c = IpData().curr_y_c()->MakeNew();
        SmartPtr<Vector> y_d = IpData().curr_y_d()->MakeNew();
        bool retval = eq_mult_calculator_->CalculateMultipliers(*y_c, *y_d);
        if (!retval) {
          y_c->Set(0.0);
          y_d->Set(0.0);
        }
        else {
          Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
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

  void MinC_1NrmRestorationPhase::ComputeBoundMultiplierStep(Vector& delta_z,
							     const Vector& curr_z,
							     const Vector& curr_slack,
							     const Vector& trial_slack)
  {
    Number mu = IpData().curr_mu();

    delta_z.Copy(curr_slack);
    delta_z.Axpy(-1., trial_slack);
    delta_z.ElementWiseMultiply(curr_z);
    delta_z.AddScalar(mu);
    delta_z.ElementWiseDivide(curr_slack);
    delta_z.Axpy(-1., curr_z);
  }

} // namespace Ipopt
