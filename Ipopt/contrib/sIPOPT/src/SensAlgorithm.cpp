// Copyright 2009 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-16

#include "SensAlgorithm.hpp"
#include "SensUtils.hpp"


namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 1;
#endif

  SensAlgorithm::SensAlgorithm(std::vector< SmartPtr<SchurDriver> >& driver_vec,
			       SmartPtr<SensitivityStepCalculator> sens_step_calc,
			       SmartPtr<Measurement> measurement,
			       Index n_sens_steps)
    :
    driver_vec_(driver_vec),
    sens_step_calc_(sens_step_calc),
    measurement_(measurement),
    n_sens_steps_(n_sens_steps) // why doesn't he get this from the options?
  {
    DBG_START_METH("SensAlgorithm::SensAlgorithm", dbg_verbosity);

    DBG_ASSERT(n_sens_steps<=driver_vec.size());
  }

  SensAlgorithm::~SensAlgorithm()
  {
    DBG_START_METH("SensAlgorithm::~SensAlgorithm", dbg_verbosity);
  }

  bool SensAlgorithm::InitializeImpl(const OptionsList& options,
				     const std::string& prefix)
  {
    return true;
  }

  /** Main loop: Wait for new measurement, Get new step, maybe deal with
   *  bounds,  see to it that everything happens in the required
   *  timeframe. */
  SensAlgorithmExitStatus SensAlgorithm::Run()
  {
    DBG_START_METH("SensAlgorithm::Run", dbg_verbosity);

    SensAlgorithmExitStatus retval = SOLVE_SUCCESS;

    /* Loop through all steps */
    SmartPtr<IteratesVector> sol = IpData().curr()->MakeNewIteratesVector();
    SmartPtr<DenseVector> delta_u;
    SmartPtr<const Vector> unscaled_x;
    SmartPtr<IteratesVector> trialcopy;
    for (Index step_i=0; step_i<n_sens_steps_; ++step_i) {
      sens_step_calc_->SetSchurDriver(driver_vec_[step_i]);
      delta_u = measurement_->GetMeasurement(step_i+1);
      delta_u->Print(Jnlst(),J_VECTOR,J_USER1,"delta_u");
      sens_step_calc_->Step(*delta_u, *sol);
      SmartPtr<IteratesVector> saved_sol = sol->MakeNewIteratesVectorCopy();
      saved_sol->Print(Jnlst(),J_VECTOR,J_USER1,"sol_vec");
      // unscale solution...
      unscaled_x = IpNLP().NLP_scaling()->unapply_vector_scaling_x(saved_sol->x());
      DBG_ASSERT(IsValid(unscaled_x));
      saved_sol->Set_x(*unscaled_x);
      unscaled_x = NULL;
      measurement_->SetSolution(step_i+1, saved_sol);

      //trialcopy = sol->MakeNewIteratesVectorCopy();
      //IpData().set_trial(trialcopy);

      /*// compute rhs (KKT evalutaion)
	SmartPtr<IteratesVector> rhs_err = IpData().curr()->MakeNewIteratesVectorCopy();
	rhs_err->Set_x(*IpCq().trial_grad_lag_x());
	rhs_err->Set_s(*IpCq().trial_grad_lag_s());
	rhs_err->Set_y_c(*IpCq().trial_c());
	rhs_err->Set_y_d(*IpCq().trial_d_minus_s());
	rhs_err->Set_z_L(*IpCq().trial_compl_x_L());
	rhs_err->Set_z_U(*IpCq().trial_compl_x_U());
	rhs_err->Set_v_L(*IpCq().trial_compl_s_L());
	rhs_err->Set_v_U(*IpCq().trial_compl_s_U());

	SmartPtr<IteratesVector> KKT_slacks;
	KKT_slacks = rhs_err->MakeNewIteratesVectorCopy();
	driver_vec_[step_i]->data_A()->TransMultiply(*delta_u, *KKT_slacks);
	KKT_slacks->Print(Jnlst(),J_VECTOR,J_USER1,"KKT_slacks");
	KKT_slacks->Axpy(1.0,*rhs_err);

	KKT_slacks->Print(Jnlst(),J_VECTOR,J_USER1,"error");

	printf("***********************************\n"
	"Running sIPOPT my-formula\n"
	"value of objective function:  %23.16e\n"
	"Nrm2 of KKT residual:         %23.16e\n"
	"Constraint violation:         %23.16e\n"
	"***********************************\n",IpCq().unscaled_trial_f(), KKT_slacks->Nrm2(), KKT_slacks->y_c()->Nrm2());
      */
      //SmartPtr<IteratesVector> trialcopyvector = IpData().curr()->MakeNewIteratesVectorCopy();
      //IpData().set_trial(trialcopyvector);
      //IpData().curr()->x()->Print(Jnlst(),J_VECTOR,J_USER1,"curr");
    }

    return retval;
  }

}
