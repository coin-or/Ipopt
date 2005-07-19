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

  DefineIpoptType(MinC_1NrmRestorationPhase);

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

  void MinC_1NrmRestorationPhase::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddLowerBoundedNumberOption(
      "bound_mult_reset_threshold",
      "Threshold for resetting bound multipliers after restoration phase.",
      0.0, false,
      1e3,
      "After returning from the restoration phase, the bound multipliers are "
      "updated with a Newton step for complementarity.  Here, the "
      "change in the primal variables during the entire restoration "
      "phase is taken to be the corresponding primal Newton step. "
      "However, if after the update the largest bound multiplier "
      "exceeds the threshold specified by this option, the multipliers "
      "are all reset to 1.");
  }

  bool MinC_1NrmRestorationPhase::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    // keep a copy of these options to use when setting up the
    // restoration phase
    resto_options_ = new OptionsList(options);

    options.GetNumericValue("constr_mult_init_max",
                            constr_mult_init_max_,
                            prefix);
    options.GetNumericValue("bound_mult_reset_threshold",
                            bound_mult_reset_threshold_,
                            prefix);
    options.GetBoolValue("expect_infeasible_problem",
                         expect_infeasible_problem_,
                         prefix);

    // ToDo take care of this somewhere else?  avoid that the
    // restoration phase is trigged by user option in first iteration
    // of the restoration phase
    resto_options_->SetValue("resto.start_with_resto", "no");

    count_restorations_ = 0;

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

    // Increase counter for restoration phase calls
    count_restorations_++;
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "Starting Restoration Phase for the %d. time\n",
                   count_restorations_);

    DBG_ASSERT(IpCq().curr_constraint_violation()>0.);

    // Create the restoration phase NLP etc objects
    SmartPtr<IpoptData> resto_ip_data = new IpoptData();
    SmartPtr<IpoptNLP> resto_ip_nlp =
      new RestoIpoptNLP(IpNLP(), IpData(), IpCq());
    SmartPtr<IpoptCalculatedQuantities> resto_ip_cq =
      new IpoptCalculatedQuantities(resto_ip_nlp, resto_ip_data);

    // Decide if we want to use the original option or want to make
    // some changes
    SmartPtr<OptionsList> actual_resto_options = resto_options_;
    if (expect_infeasible_problem_ && count_restorations_==1) {
      actual_resto_options = new OptionsList(*resto_options_);
      // Ask for significant reduction of infeasibility, in the hope
      // that we do not return from the restoration phase is the
      // problem is infeasible
      actual_resto_options->SetNumericValue("resto.kappa_resto", 1e-3);
    }

    // Initialize the restoration phase algorithm
    resto_alg_->Initialize(Jnlst(), *resto_ip_nlp, *resto_ip_data,
                           *resto_ip_cq, *actual_resto_options, "resto.");

    // Set iteration counter and info field for the restoration phase
    resto_ip_data->Set_iter_count(IpData().iter_count()+1);
    resto_ip_data->Set_info_regu_x(IpData().info_regu_x());
    resto_ip_data->Set_info_alpha_primal(IpData().info_alpha_primal());
    resto_ip_data->Set_info_alpha_primal_char(IpData().info_alpha_primal_char());
    resto_ip_data->Set_info_alpha_dual(IpData().info_alpha_dual());
    resto_ip_data->Set_info_ls_count(IpData().info_ls_count());

    // Call the optimization algorithm to solve the restoration phase
    // problem
    SolverReturn resto_status	= resto_alg_->Optimize();

    int retval=-1;

    if (resto_status == SUCCESS) {
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
                            "x", *resto_ip_data->curr()->x());
        Jnlst().PrintVector(J_VECTOR, J_LINE_SEARCH,
                            "y_c", *resto_ip_data->curr()->y_c());
        Jnlst().PrintVector(J_VECTOR, J_LINE_SEARCH,
                            "y_d", *resto_ip_data->curr()->y_d());
        Jnlst().PrintVector(J_VECTOR, J_LINE_SEARCH,
                            "z_L", *resto_ip_data->curr()->z_L());
        Jnlst().PrintVector(J_VECTOR, J_LINE_SEARCH,
                            "z_U", *resto_ip_data->curr()->z_U());
        Jnlst().PrintVector(J_VECTOR, J_LINE_SEARCH,
                            "v_L", *resto_ip_data->curr()->v_L());
        Jnlst().PrintVector(J_VECTOR, J_LINE_SEARCH,
                            "v_U", *resto_ip_data->curr()->v_U());
      }

      retval = 0;
    }
    else if (resto_status == STOP_AT_TINY_STEP ||
             resto_status == STOP_AT_ACCEPTABLE_POINT) {
      Number orig_primal_inf =
        IpCq().curr_primal_infeasibility(NORM_MAX);
      // ToDo make the factor in following line an option
      if (orig_primal_inf <= 1e2*IpData().tol()) {
        THROW_EXCEPTION(RESTORATION_FAILED,
                        "Restoration phase converged to a point with small primal infeasibility");
      }
      else {
        THROW_EXCEPTION(LOCALLY_INFEASIBLE,
                        "Restoration phase converged to a point of local infeasibility");
      }
    }
    else if (resto_status == MAXITER_EXCEEDED) {
      THROW_EXCEPTION(IpoptException, "Maximal number of iterations exceeded in restoration phase.");
      retval = 1;
    }
    else if (resto_status == LOCAL_INFEASIBILITY) {
      DBG_ASSERT(false && "Restoration of the Restoration phase should never be locally infeasible!");
      THROW_EXCEPTION(RESTORATION_FAILED, "Restoration phase itself converged to a point of local "
                      "infeasibility - This is an internal error and should not happen");
    }
    else if (resto_status == RESTORATION_FAILURE) {
      THROW_EXCEPTION(RESTORATION_FAILED, "Restoration phase in the restoration phase failed.");
    }
    else {
      Jnlst().Printf(J_ERROR, J_MAIN, "Sorry, things failed ?!?!\n");
      retval = 1;
    }

    if (retval == 0) {
      // Copy the results into the trial fields;. They will be
      // accepted later in the full algorithm
      SmartPtr<const CompoundVector> cx =
        dynamic_cast<const CompoundVector*>(GetRawPtr(resto_ip_data->curr()->x()));
      DBG_ASSERT(IsValid(cx));
      SmartPtr<IteratesVector> trial = IpData().trial()->MakeNewContainer();
      trial->Set_primal(*cx->GetComp(0), *resto_ip_data->curr()->s());
      IpData().set_trial(trial);

      // Update the bound multiplers, pretending that the entire
      // progress in x and s in the restoration phase has been one
      // [rimal-dual Newton step (and therefore the result of solving
      // an augmented system)
      SmartPtr<IteratesVector> delta = IpData().curr()->MakeNewIteratesVector(true);
      delta->Set(0.0);
      ComputeBoundMultiplierStep(*delta->z_L_NonConst(), *IpData().curr()->z_L(),
                                 *IpCq().curr_slack_x_L(),
                                 *IpCq().trial_slack_x_L());
      ComputeBoundMultiplierStep(*delta->z_U_NonConst(), *IpData().curr()->z_U(),
                                 *IpCq().curr_slack_x_U(),
                                 *IpCq().trial_slack_x_U());
      ComputeBoundMultiplierStep(*delta->v_L_NonConst(), *IpData().curr()->v_L(),
                                 *IpCq().curr_slack_s_L(),
                                 *IpCq().trial_slack_s_L());
      ComputeBoundMultiplierStep(*delta->v_U_NonConst(), *IpData().curr()->v_U(),
                                 *IpCq().curr_slack_s_U(),
                                 *IpCq().trial_slack_s_U());

      DBG_PRINT_VECTOR(1, "delta_z_L", *delta->z_L());
      DBG_PRINT_VECTOR(1, "delta_z_U", *delta->z_U());
      DBG_PRINT_VECTOR(1, "delta_v_L", *delta->v_L());
      DBG_PRINT_VECTOR(1, "delta_v_U", *delta->v_U());

      Number alpha_dual = IpCq().dual_frac_to_the_bound(IpData().curr_tau(),
                          *delta->z_L_NonConst(),
                          *delta->z_U_NonConst(),
                          *delta->v_L_NonConst(),
                          *delta->v_U_NonConst());
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "Step size for bound multipliers: %8.2e\n", alpha_dual);

      IpData().SetTrialBoundMultipliersFromStep(alpha_dual, *delta->z_L(), *delta->z_U(), *delta->v_L(), *delta->v_U() );

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
      Number bound_mult_max = Max(IpData().trial()->z_L()->Amax(),
                                  IpData().trial()->z_U()->Amax(),
                                  IpData().trial()->v_L()->Amax(),
                                  IpData().trial()->v_U()->Amax());
      trial = IpData().trial()->MakeNewContainer();
      if (bound_mult_max > bound_mult_reset_threshold_) {
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                       "Bound multipliers after restoration phase too large (max=%8.2e). Set all to 1.\n",
                       bound_mult_max);
        trial->create_new_z_L();
        trial->create_new_z_U();
        trial->create_new_v_L();
        trial->create_new_v_U();
        trial->z_L_NonConst()->Set(1.0);
        trial->z_U_NonConst()->Set(1.0);
        trial->v_L_NonConst()->Set(1.0);
        trial->v_U_NonConst()->Set(1.0);
        IpData().set_trial(trial);

      }
      // Recompute the equality constraint multipliers as least square estimate
      trial = IpData().trial()->MakeNewContainer();
      if (IsValid(eq_mult_calculator_) && constr_mult_init_max_>0.) {
        // First move all the trial data into the current fields, since
        // those values are needed to compute the initial values for
        // the multipliers
        IpData().CopyTrialToCurrent();
        SmartPtr<Vector> y_c = IpData().curr()->y_c()->MakeNew();
        SmartPtr<Vector> y_d = IpData().curr()->y_d()->MakeNew();
        bool retval = eq_mult_calculator_->CalculateMultipliers(*y_c, *y_d);
        if (!retval) {
          y_c->Set(0.0);
          y_d->Set(0.0);
        }
        else {
          Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                         "Least square estimates max(y_c) = %e, max(y_d) = %e\n",
                         y_c->Amax(), y_d->Amax());
          Number y_init_nrm = Max(y_c->Amax(), y_d->Amax());
          if (!retval || y_init_nrm > constr_mult_init_max_) {
            y_c->Set(0.0);
            y_d->Set(0.0);
          }
        }
        trial->Set_eq_mult(*y_c, *y_d);
      }
      else {
        SmartPtr<Vector> y_c = IpData().curr()->y_c()->MakeNew();
        SmartPtr<Vector> y_d = IpData().curr()->y_d()->MakeNew();
        y_c->Set(0.0);
        y_d->Set(0.0);
        trial->Set_eq_mult(*y_c, *y_d);
      }
      IpData().set_trial(trial);

      DBG_PRINT_VECTOR(2, "y_c", *IpData().curr()->y_c());
      DBG_PRINT_VECTOR(2, "y_d", *IpData().curr()->y_d());

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
