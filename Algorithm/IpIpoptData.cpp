// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpIpoptData.hpp"
#include "IpIpoptNLP.hpp"

namespace Ipopt
{

  DBG_SET_VERBOSITY(0);

  IpoptData::IpoptData()
      :
      iter_count_(0),
      curr_mu_(-1.),
      mu_initialized_(false),
      curr_tau_(-1.),
      tau_initialized_(false),
      initialize_called_(false),
      have_prototypes_(false),
      free_mu_mode_(false),

      info_alpha_primal_(0.),
      info_alpha_primal_char_(' '),
      info_alpha_dual_(0.),
      info_regu_x_(0.),
      info_ls_count_(0),
      info_skip_output_(false)
  {}

  IpoptData::~IpoptData()
  {}

  bool IpoptData::Initialize(const Journalist& jnlst,
                             const OptionsList& options,
                             const std::string& prefix)
  {
    Number value;

    if (options.GetNumericValue("epsilon_tol", value, prefix)) {
      ASSERT_EXCEPTION(value > 0.0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"epsilon_tol\": This value must be larger than 0.");
      epsilon_tol_ = value;
    }
    else {
      epsilon_tol_ = 1e-8;
    }

    iter_count_=0;

    have_prototypes_ = false;
    tau_initialized_ = false;
    have_deltas_ = false;

    initialize_called_ = true;
    return true;
  }

  bool IpoptData::InitializeDataStructures(IpoptNLP& ip_nlp,
      bool want_x,
      bool want_y_c,
      bool want_y_d,
      bool want_z_L,
      bool want_z_U,
      bool want_v_L,
      bool want_v_U)
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
                                  new_v_L, want_v_L,
                                  new_v_U, want_v_U);
    if (!retValue) {
      return false;
    }

    new_s = new_y_d->MakeNew(); // same dimension as d
    //new_s_->Set(0.0);

    curr_x_ = ConstPtr(new_x);
    curr_s_ = ConstPtr(new_s);
    curr_y_c_ = ConstPtr(new_y_c);
    curr_y_d_ = ConstPtr(new_y_d);
    curr_z_L_ = ConstPtr(new_z_L);
    curr_z_U_ = ConstPtr(new_z_U);
    curr_v_L_ = ConstPtr(new_v_L);
    curr_v_U_ = ConstPtr(new_v_U);

    // create structures for the deltas (uninitialized)
    delta_x_ = curr_x_->MakeNew();
    delta_s_ = curr_s_->MakeNew();
    delta_y_c_ = curr_y_c_->MakeNew();
    delta_y_d_ = curr_y_d_->MakeNew();
    delta_z_L_ = curr_z_L_->MakeNew();
    delta_z_U_ = curr_z_U_->MakeNew();
    delta_v_L_ = curr_v_L_->MakeNew();
    delta_v_U_ = curr_v_U_->MakeNew();

    // Initialize the steps to 0 (so that we can access them for
    // example in IterationOutput before the first iteration)
    delta_x_->Set(0.0);
    delta_s_->Set(0.0);
    delta_y_c_->Set(0.0);
    delta_y_d_->Set(0.0);
    delta_z_L_->Set(0.0);
    delta_z_U_->Set(0.0);
    delta_v_L_->Set(0.0);
    delta_v_U_->Set(0.0);

    have_prototypes_ = true;
    have_deltas_ = false;

    return true;
  }

  void IpoptData::SetTrialPrimalVariables(const Vector& x, const Vector& s)
  {
    DBG_START_METH("IpoptData::SetTrialPrimalVariables", 0);
    DBG_ASSERT(have_prototypes_);

    SmartPtr<Vector> newvec = curr_x_->MakeNew();
    newvec->Copy(x);
    trial_x_ = ConstPtr(newvec);

    newvec = curr_s_->MakeNew();
    newvec->Copy(s);
    trial_s_ = ConstPtr(newvec);
  }

  void IpoptData::SetTrialXVariables(const Vector& x)
  {
    DBG_START_METH("IpoptData::SetTrialXVariables", 0);
    DBG_ASSERT(have_prototypes_);

    SmartPtr<Vector> newvec = curr_x_->MakeNew();
    newvec->Copy(x);
    trial_x_ = ConstPtr(newvec);
  }

  void IpoptData::SetTrialSVariables(const Vector& s)
  {
    DBG_START_METH("IpoptData::SetTrialSVariables", 0);
    DBG_ASSERT(have_prototypes_);

    SmartPtr<Vector> newvec = curr_s_->MakeNew();
    newvec->Copy(s);
    trial_s_ = ConstPtr(newvec);
  }

  void IpoptData::SetTrialEqMultipliers(const Vector& y_c, const Vector& y_d)
  {
    DBG_ASSERT(have_prototypes_);

    SmartPtr<Vector> newvec = curr_y_c_->MakeNew();
    newvec->Copy(y_c);
    trial_y_c_=ConstPtr(newvec);

    newvec = curr_y_d_->MakeNew();
    newvec->Copy(y_d);
    trial_y_d_=ConstPtr(newvec);
  }

  void IpoptData::SetTrialBoundMultipliers(const Vector& z_L,
      const Vector& z_U,
      const Vector& v_L,
      const Vector& v_U)
  {
    DBG_ASSERT(have_prototypes_);

    SmartPtr<Vector> newvec = curr_z_L_->MakeNew();
    newvec->Copy(z_L);
    trial_z_L_=ConstPtr(newvec);

    newvec = curr_z_U_->MakeNew();
    newvec->Copy(z_U);
    trial_z_U_=ConstPtr(newvec);

    newvec = curr_v_L_->MakeNew();
    newvec->Copy(v_L);
    trial_v_L_=ConstPtr(newvec);

    newvec = curr_v_U_->MakeNew();
    newvec->Copy(v_U);
    trial_v_U_=ConstPtr(newvec);
  }

  void IpoptData::SetTrialPrimalVariablesFromStep(Number alpha,
      const Vector& delta_x,
      const Vector& delta_s)
  {
    DBG_ASSERT(have_prototypes_);

    SmartPtr<Vector> newvec = curr_x_->MakeNew();
    newvec->AddTwoVectors(1., *curr_x_, alpha, delta_x, 0.);
    /* DELE
    newvec->Copy(*curr_x_);
    newvec->Axpy(alpha, delta_x);
    */
    trial_x_ = ConstPtr(newvec);

    newvec = curr_s_->MakeNew();
    newvec->Copy(*curr_s_);
    newvec->Axpy(alpha, delta_s);
    trial_s_ = ConstPtr(newvec);
  }

  void IpoptData::SetTrialEqMultipilersFromStep(Number alpha,
      const Vector& delta_y_c,
      const Vector& delta_y_d)
  {
    DBG_ASSERT(have_prototypes_);

    SmartPtr<Vector> newvec = curr_y_c_->MakeNew();
    newvec->Copy(*curr_y_c_);
    newvec->Axpy(alpha, delta_y_c);
    trial_y_c_ = ConstPtr(newvec);

    newvec = curr_y_d_->MakeNew();
    newvec->Copy(*curr_y_d_);
    newvec->Axpy(alpha, delta_y_d);
    trial_y_d_ = ConstPtr(newvec);
  }

  void IpoptData::SetTrialBoundMultipliersFromStep(Number alpha,
      const Vector& delta_z_L,
      const Vector& delta_z_U,
      const Vector& delta_v_L,
      const Vector& delta_v_U)
  {
    DBG_ASSERT(have_prototypes_);

    SmartPtr<Vector> newvec = curr_z_L_->MakeNew();
    newvec->Copy(*curr_z_L_);
    newvec->Axpy(alpha, delta_z_L);
    trial_z_L_ = ConstPtr(newvec);

    newvec = curr_z_U_->MakeNew();
    newvec->Copy(*curr_z_U_);
    newvec->Axpy(alpha, delta_z_U);
    trial_z_U_ = ConstPtr(newvec);

    newvec = curr_v_L_->MakeNew();
    newvec->Copy(*curr_v_L_);
    newvec->Axpy(alpha, delta_v_L);
    trial_v_L_ = ConstPtr(newvec);

    newvec = curr_v_U_->MakeNew();
    newvec->Copy(*curr_v_U_);
    newvec->Axpy(alpha, delta_v_U);
    trial_v_U_ = ConstPtr(newvec);
  }

  void IpoptData::SetTrialPrimalVariablesFromPtr(
    const SmartPtr<const Vector>& xptr,
    const SmartPtr<const Vector>& sptr)
  {
    DBG_ASSERT(have_prototypes_);

    trial_x_ = xptr;
    trial_s_ = sptr;
  }

  void IpoptData::SetTrialConstraintMultipliersFromPtr(
    const SmartPtr<const Vector>& y_cptr,
    const SmartPtr<const Vector>& y_dptr)
  {
    DBG_ASSERT(have_prototypes_);

    trial_y_c_ = y_cptr;
    trial_y_d_ = y_dptr;
  }

  void IpoptData::SetTrialBoundMultipliersFromPtr(
    const SmartPtr<const Vector>& z_Lptr,
    const SmartPtr<const Vector>& z_Uptr,
    const SmartPtr<const Vector>& v_Lptr,
    const SmartPtr<const Vector>& v_Uptr)
  {
    DBG_ASSERT(have_prototypes_);

    trial_z_L_ = z_Lptr;
    trial_z_U_ = z_Uptr;
    trial_v_L_ = v_Lptr;
    trial_v_U_ = v_Uptr;
  }

  void IpoptData::CopyTrialToCurrent()
  {
    curr_x_ = trial_x_;
    curr_s_ = trial_s_;
    curr_y_c_ = trial_y_c_;
    curr_y_d_ = trial_y_d_;
    curr_z_L_ = trial_z_L_;
    curr_z_U_ = trial_z_U_;
    curr_v_L_ = trial_v_L_;
    curr_v_U_ = trial_v_U_;
  }

  void IpoptData::AcceptTrialPoint()
  {
    DBG_ASSERT(IsValid(trial_x_));
    DBG_ASSERT(IsValid(trial_s_));
    DBG_ASSERT(IsValid(trial_y_c_));
    DBG_ASSERT(IsValid(trial_y_d_));
    DBG_ASSERT(IsValid(trial_z_L_));
    DBG_ASSERT(IsValid(trial_z_U_));
    DBG_ASSERT(IsValid(trial_v_L_));
    DBG_ASSERT(IsValid(trial_v_U_));

    CopyTrialToCurrent();

    // Set trial pointers to Null (frees memory unless someone else is
    // still referring to it, and it makes sure that indeed all trial
    // values are set before a new trial point is accepted)
    trial_x_ = NULL;
    trial_s_ = NULL;
    trial_y_c_ = NULL;
    trial_y_d_ = NULL;
    trial_z_L_ = NULL;
    trial_z_U_ = NULL;
    trial_v_L_ = NULL;
    trial_v_U_ = NULL;

    // ToDo that should be handled better
    have_deltas_ = false;
  }

} // namespace Ipopt
