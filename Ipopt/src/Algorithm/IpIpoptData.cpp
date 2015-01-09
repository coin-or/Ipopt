// Copyright (C) 2004, 2010 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpIpoptData.hpp"
#include "IpIpoptNLP.hpp"

namespace Ipopt
{

  IpoptData::IpoptData(SmartPtr<IpoptAdditionalData> add_data /*= NULL*/,
                       Number cpu_time_start /*= -1.*/)
      :
      cpu_time_start_(cpu_time_start),
      add_data_(add_data)
  {}

  IpoptData::~IpoptData()
  {}

  void IpoptData::RegisterOptions(const SmartPtr<RegisteredOptions>& roptions)
  {
    roptions->SetRegisteringCategory("Convergence");
    roptions->AddLowerBoundedNumberOption(
      "tol",
      "Desired convergence tolerance (relative).",
      0.0, true,  1e-8,
      "Determines the convergence tolerance for the algorithm.  The "
      "algorithm terminates successfully, if the (scaled) NLP error "
      "becomes smaller than this value, and if the (absolute) criteria "
      "according to \"dual_inf_tol\", \"constr_viol_tol\", and "
      "\"compl_inf_tol\" are met.  (This is epsilon_tol in Eqn. (6) in "
      "implementation paper).  See also \"acceptable_tol\" as a second "
      "termination criterion.  Note, some other algorithmic features also use "
      "this quantity to determine thresholds etc.");
  }

  bool IpoptData::Initialize(const Journalist& jnlst,
                             const OptionsList& options,
                             const std::string& prefix)
  {
#if 0
    // I (AW) took the following heuristic out again, since it seemed
    // that the restoration phase tolerance became too tight by
    // default.  I originally probably put it in to avoid that a claim
    // of infeasibility is made prematurely...  let's see if someone
    // starts screaming...
    if (prefix=="resto.") {
      // The default for the restoration phase is 1e-2 time the value
      // for the regular algorithm
      if (!options.GetNumericValue("resto.tol", tol_, "")) {
        options.GetNumericValue("tol", tol_, prefix);
        tol_ *= 1e-2;
      }
    }
    else {
      options.GetNumericValue("tol", tol_, prefix);
    }
#else
    options.GetNumericValue("tol", tol_, prefix);
#endif

    iter_count_ = 0;
    curr_mu_ = -1.;
    mu_initialized_ = false;
    curr_tau_ = -1.;
    tau_initialized_ = false;
    have_prototypes_ = false;
    have_deltas_ = false;
    have_affine_deltas_ = false;

    free_mu_mode_ = false;
    tiny_step_flag_ = false;

    info_ls_count_ = 0;
    ResetInfo();
    info_last_output_ = -1.;
    info_iters_since_header_ = 1000; // need to be larger 10

    initialize_called_ = true;

    if (cpu_time_start_==-1.) {
      cpu_time_start_ = CpuTime();
    }

    bool retval = true;

    if (IsValid(add_data_)) {
      retval = add_data_->Initialize(jnlst, options, prefix);
    }

    return retval;
  }

  bool IpoptData::InitializeDataStructures(IpoptNLP& ip_nlp,
      bool want_x,
      bool want_y_c,
      bool want_y_d,
      bool want_z_L,
      bool want_z_U)
  {
    DBG_ASSERT(initialize_called_);
    /*
     * Allocate space for all the required linear algebra 
     * structures
     */

    SmartPtr<Vector> new_x;
    SmartPtr<Vector> new_s;
    SmartPtr<Vector> new_y_c;
    SmartPtr<Vector> new_y_d;
    SmartPtr<Vector> new_z_L;
    SmartPtr<Vector> new_z_U;
    SmartPtr<Vector> new_v_L;
    SmartPtr<Vector> new_v_U;

    // Get the required linear algebra structures from the model
    bool retValue
    = ip_nlp.InitializeStructures(new_x, want_x,
                                  new_y_c, want_y_c,
                                  new_y_d, want_y_d,
                                  new_z_L, want_z_L,
                                  new_z_U, want_z_U,
                                  new_v_L, new_v_U);
    if (!retValue) {
      return false;
    }

    new_s = new_y_d->MakeNew(); // same dimension as d

    iterates_space_ = new IteratesVectorSpace(*(new_x->OwnerSpace()), *(new_s->OwnerSpace()),
                      *(new_y_c->OwnerSpace()), *(new_y_d->OwnerSpace()),
                      *(new_z_L->OwnerSpace()), *(new_z_U->OwnerSpace()),
                      *(new_v_L->OwnerSpace()), *(new_v_U->OwnerSpace())
                                             );

    curr_ = iterates_space_->MakeNewIteratesVector(*new_x,
            *new_s,
            *new_y_c,
            *new_y_d,
            *new_z_L,
            *new_z_U,
            *new_v_L,
            *new_v_U);
#if COIN_IPOPT_CHECKLEVEL > 0

    debug_curr_tag_ = curr_->GetTag();
    debug_curr_tag_sum_ = curr_->GetTagSum();
    debug_trial_tag_ = 0;
    debug_trial_tag_sum_ = 0;
    debug_delta_tag_ = 0;
    debug_delta_tag_sum_ = 0;
    debug_delta_aff_tag_ = 0;
    debug_delta_aff_tag_sum_ = 0;
#endif

    trial_ = NULL;

    // Set the pointers for storing steps to NULL
    delta_ = NULL;
    delta_aff_ = NULL;

    have_prototypes_ = true;
    have_deltas_ = false;
    have_affine_deltas_ = false;

    bool retval = true;
    if (IsValid(add_data_)) {
      retval = add_data_->InitializeDataStructures();
    }

    return retval;
  }

  void IpoptData::SetTrialPrimalVariablesFromStep(Number alpha,
      const Vector& delta_x,
      const Vector& delta_s)
  {
    DBG_ASSERT(have_prototypes_);

    if (IsNull(trial_)) {
      trial_ = iterates_space_->MakeNewIteratesVector(false);
    }

    SmartPtr<IteratesVector> newvec = trial_->MakeNewContainer();
    newvec->create_new_x();
    newvec->x_NonConst()->AddTwoVectors(1., *curr_->x(), alpha, delta_x, 0.);

    newvec->create_new_s();
    newvec->s_NonConst()->AddTwoVectors(1., *curr_->s(), alpha, delta_s, 0.);

    set_trial(newvec);
  }

  void IpoptData::SetTrialEqMultipliersFromStep(Number alpha,
      const Vector& delta_y_c,
      const Vector& delta_y_d)
  {
    DBG_ASSERT(have_prototypes_);

    SmartPtr<IteratesVector> newvec = trial()->MakeNewContainer();
    newvec->create_new_y_c();
    newvec->y_c_NonConst()->AddTwoVectors(1., *curr()->y_c(), alpha, delta_y_c, 0.);

    newvec->create_new_y_d();
    newvec->y_d_NonConst()->AddTwoVectors(1., *curr()->y_d(), alpha, delta_y_d, 0.);

    set_trial(newvec);
  }

  void IpoptData::SetTrialBoundMultipliersFromStep(Number alpha,
      const Vector& delta_z_L,
      const Vector& delta_z_U,
      const Vector& delta_v_L,
      const Vector& delta_v_U)
  {
    DBG_ASSERT(have_prototypes_);

    SmartPtr<IteratesVector> newvec = trial()->MakeNewContainer();
    newvec->create_new_z_L();
    newvec->z_L_NonConst()->AddTwoVectors(1., *curr()->z_L(), alpha, delta_z_L, 0.);

    newvec->create_new_z_U();
    newvec->z_U_NonConst()->AddTwoVectors(1., *curr()->z_U(), alpha, delta_z_U, 0.);

    newvec->create_new_v_L();
    newvec->v_L_NonConst()->AddTwoVectors(1., *curr()->v_L(), alpha, delta_v_L, 0.);

    newvec->create_new_v_U();
    newvec->v_U_NonConst()->AddTwoVectors(1., *curr()->v_U(), alpha, delta_v_U, 0.);

    set_trial(newvec);
  }

  void IpoptData::AcceptTrialPoint()
  {
    DBG_ASSERT(IsValid(trial_));
    DBG_ASSERT(IsValid(trial_->x()));
    DBG_ASSERT(IsValid(trial_->s()));
    DBG_ASSERT(IsValid(trial_->y_c()));
    DBG_ASSERT(IsValid(trial_->y_d()));
    DBG_ASSERT(IsValid(trial_->z_L()));
    DBG_ASSERT(IsValid(trial_->z_U()));
    DBG_ASSERT(IsValid(trial_->v_L()));
    DBG_ASSERT(IsValid(trial_->v_U()));

    CopyTrialToCurrent();

    // Set trial pointers to Null (frees memory unless someone else is
    // still referring to it, and it makes sure that indeed all trial
    // values are set before a new trial point is accepted)
    trial_ = NULL;

    // Free the memory for the affine-scaling step
    delta_aff_ = NULL;

    have_deltas_ = false;
    have_affine_deltas_ = false;

    if (IsValid(add_data_)) {
      add_data_->AcceptTrialPoint();
    }
  }

} // namespace Ipopt
