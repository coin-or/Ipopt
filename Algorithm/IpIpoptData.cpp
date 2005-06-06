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
      tiny_step_flag_(false),

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

    if (options.GetNumericValue("tol", value, prefix)) {
      ASSERT_EXCEPTION(value > 0.0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"tol\": This value must be larger than 0.");
      tol_ = value;
    }
    else {
      tol_ = 1e-8;
    }

    if (options.GetNumericValue("dual_inf_tol", value, prefix)) {
      ASSERT_EXCEPTION(value > 0., OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"dual_inf_tol\": This value must be larger than 0.");
      dual_inf_tol_ = value;
    }
    else {
      dual_inf_tol_ = 1e-2;
    }

    if (options.GetNumericValue("primal_inf_tol", value, prefix)) {
      ASSERT_EXCEPTION(value > 0., OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"primal_inf_tol\": This value must be larger than 0.");
      primal_inf_tol_ = value;
    }
    else {
      primal_inf_tol_ = 1e-2;
    }

    if (options.GetNumericValue("compl_inf_tol", value, prefix)) {
      ASSERT_EXCEPTION(value > 0., OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"compl_inf_tol\": This value must be larger than 0.");
      compl_inf_tol_ = value;
    }
    else {
      compl_inf_tol_ = 1e-2;
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

    trial_ = NULL;

    // Set the pointers for storing steps to NULL
    delta_x_ = NULL;
    delta_s_ = NULL;
    delta_y_c_ = NULL;
    delta_y_d_ = NULL;
    delta_z_L_ = NULL;
    delta_z_U_ = NULL;
    delta_v_L_ = NULL;
    delta_v_U_ = NULL;

    // Set the pointers for storing steps to NULL
    delta_aff_x_ = NULL;
    delta_aff_s_ = NULL;
    delta_aff_y_c_ = NULL;
    delta_aff_y_d_ = NULL;
    delta_aff_z_L_ = NULL;
    delta_aff_z_U_ = NULL;
    delta_aff_v_L_ = NULL;
    delta_aff_v_U_ = NULL;

    have_prototypes_ = true;
    have_deltas_ = false;

    return true;
  }

  SmartPtr<const IteratesVector> IpoptData::curr() const
  {
    return curr_;
//     // create a new container
//     SmartPtr<IteratesVector> iterates = iterates_space_->MakeNewIteratesVector(false);
//     if (IsValid(curr_x_)) { iterates->Set_x(*curr_x_); }
//     if (IsValid(curr_s_)) { iterates->Set_s(*curr_s_); }
//     if (IsValid(curr_y_c_)) { iterates->Set_y_c(*curr_y_c_); }
//     if (IsValid(curr_y_d_)) { iterates->Set_y_d(*curr_y_d_); }
//     if (IsValid(curr_z_L_)) { iterates->Set_z_L(*curr_z_L_); }
//     if (IsValid(curr_z_U_)) { iterates->Set_z_U(*curr_z_U_); }
//     if (IsValid(curr_v_L_)) {iterates->Set_v_L(*curr_v_L_); }
//     if (IsValid(curr_v_U_)) { iterates->Set_v_U(*curr_v_U_); }

//     return ConstPtr(iterates);
  }

  SmartPtr<const IteratesVector> IpoptData::trial() const
  {
#ifdef IP_DEBUG
    DBG_ASSERT(IsNull(trial_) || trial_->GetTag() == debug_trial_tag_);
#endif
    return ConstPtr(trial_);

// #ifdef IP_DEBUG
//     // A debug check to make sure that the trial data was not changed behind
//     // IpoptData's back.
//     DBG_ASSERT((IsNull(trial_x_) || trial_x_->GetTag() == debug_trial_x_tag_)
// 	       && (IsNull(trial_s_) || trial_s_->GetTag() == debug_trial_s_tag_)
// 	       && (IsNull(trial_y_c_) || trial_y_c_->GetTag() == debug_trial_y_c_tag_)
// 	       && (IsNull(trial_y_d_) || trial_y_d_->GetTag() == debug_trial_y_d_tag_)
// 	       && (IsNull(trial_z_L_) || trial_z_L_->GetTag() == debug_trial_z_L_tag_)
// 	       && (IsNull(trial_z_U_) || trial_z_U_->GetTag() == debug_trial_z_U_tag_)
// 	       && (IsNull(trial_v_L_) || trial_v_L_->GetTag() == debug_trial_v_L_tag_)
// 	       && (IsNull(trial_v_U_) || trial_v_U_->GetTag() == debug_trial_v_U_tag_)
// 	       );
// #endif

//     // for now, we create a new container - when this class uses an IteratesVector
//     // internally, we wil be able to efficiently and safely return that vector.
//     SmartPtr<IteratesVector> iterates = iterates_space_->MakeNewIteratesVector(false);
//     if (IsValid(trial_x_)) { iterates->Set_x(*trial_x_); }
//     if (IsValid(trial_s_)) { iterates->Set_s(*trial_s_); }
//     if (IsValid(trial_y_c_)) { iterates->Set_y_c(*trial_y_c_); }
//     if (IsValid(trial_y_d_)) { iterates->Set_y_d(*trial_y_d_); }
//     if (IsValid(trial_z_L_)) { iterates->Set_z_L(*trial_z_L_); }
//     if (IsValid(trial_z_U_)) { iterates->Set_z_U(*trial_z_U_); }
//     if (IsValid(trial_v_L_)) { iterates->Set_v_L(*trial_v_L_); }
//     if (IsValid(trial_v_U_)) { iterates->Set_v_U(*trial_v_U_); }

//     return ConstPtr(iterates);
  }

  SmartPtr<const IteratesVector> IpoptData::delta() const
  {
    //ToDo: Add the debug tag check in here...

    // create a new container
    SmartPtr<IteratesVector> iterates = iterates_space_->MakeNewIteratesVector(false);
    if (IsValid(delta_x_)) { iterates->Set_x(*delta_x_); }
    if (IsValid(delta_s_)) { iterates->Set_s(*delta_s_); }
    if (IsValid(delta_y_c_)) { iterates->Set_y_c(*delta_y_c_); }
    if (IsValid(delta_y_d_)) { iterates->Set_y_d(*delta_y_d_); }
    if (IsValid(delta_z_L_)) { iterates->Set_z_L(*delta_z_L_); }
    if (IsValid(delta_z_U_)) { iterates->Set_z_U(*delta_z_U_); }
    if (IsValid(delta_v_L_)) { iterates->Set_v_L(*delta_v_L_); }
    if (IsValid(delta_v_U_)) { iterates->Set_v_U(*delta_v_U_); }

    return ConstPtr(iterates);
  }

  SmartPtr<const IteratesVector> IpoptData::delta_aff() const
  {
    //ToDo: Add the debug tag check in here
    
    // create a new container
    SmartPtr<IteratesVector> iterates = iterates_space_->MakeNewIteratesVector(false);
    if (IsValid(delta_aff_x_)) { iterates->Set_x(*delta_aff_x_); }
    if (IsValid(delta_aff_s_)) { iterates->Set_s(*delta_aff_s_); }
    if (IsValid(delta_aff_y_c_)) { iterates->Set_y_c(*delta_aff_y_c_); }
    if (IsValid(delta_aff_y_d_)) { iterates->Set_y_d(*delta_aff_y_d_); }
    if (IsValid(delta_aff_z_L_)) { iterates->Set_z_L(*delta_aff_z_L_); }
    if (IsValid(delta_aff_z_U_)) { iterates->Set_z_U(*delta_aff_z_U_); }
    if (IsValid(delta_aff_v_L_)) { iterates->Set_v_L(*delta_aff_v_L_); }
    if (IsValid(delta_aff_v_U_)) { iterates->Set_v_U(*delta_aff_v_U_); }

    return ConstPtr(iterates);
  }

  void IpoptData::set_trial(const SmartPtr<const IteratesVector>& trial)
  {
    trial_ = trial;
#ifdef IP_DEBUG
    if (IsValid(trial)) {
      debug_trial_tag_ = trial->GetTag();
    }
    else {
      debug_trial_tag_ = 0;
    }
#endif
//     trial_x_ = trial->x();
//     trial_s_ = trial->s();
//     trial_y_c_ = trial->y_c();
//     trial_y_d_ = trial->y_d();
//     trial_z_L_ = trial->z_L();
//     trial_z_U_ = trial->z_U();
//     trial_v_L_ = trial->v_L();
//     trial_v_U_ = trial->v_U();

// #ifdef IP_DEBUG
//     // This is to allow a debug check to make sure that data is not changed
//     // behind IpoptData's back
//     debug_trial_x_tag_ = (IsValid(trial_x_)) ? trial_x_->GetTag() : 0;
//     debug_trial_s_tag_ = (IsValid(trial_s_)) ? trial_s_->GetTag() : 0;
//     debug_trial_y_c_tag_ = (IsValid(trial_y_c_)) ? trial_y_c_->GetTag() : 0;
//     debug_trial_y_d_tag_ = (IsValid(trial_y_d_)) ? trial_y_d_->GetTag() : 0;
//     debug_trial_z_L_tag_ = (IsValid(trial_z_L_)) ? trial_z_L_->GetTag() : 0;
//     debug_trial_z_U_tag_ = (IsValid(trial_z_U_)) ? trial_z_U_->GetTag() : 0;
//     debug_trial_v_L_tag_ = (IsValid(trial_v_L_)) ? trial_v_L_->GetTag() : 0;
//     debug_trial_v_U_tag_ = (IsValid(trial_v_U_)) ? trial_v_U_->GetTag() : 0;
// #endif
  }
    
  void IpoptData::set_delta(SmartPtr<IteratesVector>& delta)
  {
    delta_x_ = delta->x();
    delta_s_ = delta->s();
    delta_y_c_ = delta->y_c();
    delta_y_d_ = delta->y_d();
    delta_z_L_ = delta->z_L();
    delta_z_U_ = delta->z_U();
    delta_v_L_ = delta->v_L();
    delta_v_U_ = delta->v_U();
  }

  void IpoptData::set_delta_aff(SmartPtr<IteratesVector>& delta_aff)
  {
    delta_aff_x_ = delta_aff->x();
    delta_aff_s_ = delta_aff->s();
    delta_aff_y_c_ = delta_aff->y_c();
    delta_aff_y_d_ = delta_aff->y_d();
    delta_aff_z_L_ = delta_aff->z_L();
    delta_aff_z_U_ = delta_aff->z_U();
    delta_aff_v_L_ = delta_aff->v_L();
    delta_aff_v_U_ = delta_aff->v_U();
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

//   void IpoptData::SetTrialFromStep(Number primal_alpha, 
// 				   Number eq_mult_alpha,
// 				   Number bound_mult_alpha,
// 				   const SmartPtr<const IteratesVector>& delta)
//   {
//     DBG_ASSERT(have_prototypes_);

//     // Once IpData stores an iterate, this method will change
//     // and use the true trial_iterate, not create a container
//     // and set...
//     SmartPtr<IteratesVector> new_trial = curr()->MakeNewContainer();
//     if (primal_alpha != 0) {
//       new_trial->create_new_x();
//       new_trial->x_NonConst()->AddTwoVectors(1.0, *curr()->x(), primal_alpha, *delta->x(), 0.0);

//       new_trial->create_new_s();
//       new_trial->s_NonConst()->AddTwoVectors(1.0, *curr()->s(), primal_alpha, *delta->s(), 0.0);
//     }

//     if (eq_mult_alpha != 0) {
//       new_trial->create_new_y_c();
//       new_trial->y_c_NonConst()->AddTwoVectors(1.0, *curr()->y_c(), eq_mult_alpha, *delta->y_c(), 0.0);
      
//       new_trial->create_new_y_d();
//       new_trial->y_d_NonConst()->AddTwoVectors(1.0, *curr()->y_d(), eq_mult_alpha, *delta->y_d(), 0.0);
//     }

//     if (bound_mult_alpha != 0) {
//       new_trial->create_new_z_L();
//       new_trial->z_L_NonConst()->AddTwoVectors(1.0, *curr()->z_L(), bound_mult_alpha, *delta->z_L(), 0.0);

//       new_trial->create_new_z_U();
//       new_trial->z_U_NonConst()->AddTwoVectors(1.0, *curr()->z_U(), bound_mult_alpha, *delta->z_U(), 0.0);

//       new_trial->create_new_v_L();
//       new_trial->v_L_NonConst()->AddTwoVectors(1.0, *curr()->v_L(), bound_mult_alpha, *delta->v_L(), 0.0);

//       new_trial->create_new_v_U();
//       new_trial->v_U_NonConst()->AddTwoVectors(1.0, *curr()->v_U(), bound_mult_alpha, *delta->v_U(), 0.0);
//     }

//     set_trial(new_trial);
//   }

//   void IpoptData::SetTrialPrimalVariablesFromStep(Number alpha,
//       const Vector& delta_x,
//       const Vector& delta_s)
//   {
//     DBG_ASSERT(have_prototypes_);

//     SmartPtr<Vector> newvec = curr_x_->MakeNew();
//     newvec->AddTwoVectors(1., *curr_x_, alpha, delta_x, 0.);
//     /* DELE
//     newvec->Copy(*curr_x_);
//     newvec->Axpy(alpha, delta_x);
//     */
//     trial_x_ = ConstPtr(newvec);

//     newvec = curr_s_->MakeNew();
//     newvec->AddTwoVectors(1., *curr_s_, alpha, delta_s, 0.);
//     /* DELE
//     newvec->Copy(*curr_s_);
//     newvec->Axpy(alpha, delta_s);
//     */
//     trial_s_ = ConstPtr(newvec);
//   }

//   void IpoptData::SetTrialEqMultipilersFromStep(Number alpha,
//       const Vector& delta_y_c,
//       const Vector& delta_y_d)
//   {
//     DBG_ASSERT(have_prototypes_);

//     SmartPtr<Vector> newvec = curr_y_c_->MakeNew();
//     newvec->AddTwoVectors(1., *curr_y_c_, alpha, delta_y_c, 0.);
//     /* DELE
//     newvec->Copy(*curr_y_c_);
//     newvec->Axpy(alpha, delta_y_c);
//     */
//     trial_y_c_ = ConstPtr(newvec);

//     newvec = curr_y_d_->MakeNew();
//     newvec->AddTwoVectors(1., *curr_y_d_, alpha, delta_y_d, 0.);
//     /* DELE
//     newvec->Copy(*curr_y_d_);
//     newvec->Axpy(alpha, delta_y_d);
//     */
//     trial_y_d_ = ConstPtr(newvec);
//   }

//   void IpoptData::SetTrialBoundMultipliersFromStep(Number alpha,
//       const Vector& delta_z_L,
//       const Vector& delta_z_U,
//       const Vector& delta_v_L,
//       const Vector& delta_v_U)
//   {
//     DBG_ASSERT(have_prototypes_);

//     SmartPtr<Vector> newvec = curr_z_L_->MakeNew();
//     newvec->AddTwoVectors(1., *curr_z_L_, alpha, delta_z_L, 0.);
//     /* DELE
//     newvec->Copy(*curr_z_L_);
//     newvec->Axpy(alpha, delta_z_L);
//     */
//     trial_z_L_ = ConstPtr(newvec);

//     newvec = curr_z_U_->MakeNew();
//     newvec->AddTwoVectors(1., *curr_z_U_, alpha, delta_z_U, 0.);
//     /* DELE
//     newvec->Copy(*curr_z_U_);
//     newvec->Axpy(alpha, delta_z_U);
//     */
//     trial_z_U_ = ConstPtr(newvec);

//     newvec = curr_v_L_->MakeNew();
//     newvec->AddTwoVectors(1., *curr_v_L_, alpha, delta_v_L, 0.);
//     /* DELE
//     newvec->Copy(*curr_v_L_);
//     newvec->Axpy(alpha, delta_v_L);
//     */
//     trial_v_L_ = ConstPtr(newvec);

//     newvec = curr_v_U_->MakeNew();
//     newvec->AddTwoVectors(1., *curr_v_U_, alpha, delta_v_U, 0.);
//     /* DELE
//     newvec->Copy(*curr_v_U_);
//     newvec->Axpy(alpha, delta_v_U);
//     */
//     trial_v_U_ = ConstPtr(newvec);
//   }

//   void IpoptData::SetTrialPrimalVariablesFromPtr(
//     const SmartPtr<const Vector>& xptr,
//     const SmartPtr<const Vector>& sptr)
//   {
//     DBG_ASSERT(have_prototypes_);

//     trial_x_ = xptr;
//     trial_s_ = sptr;
//   }

//   void IpoptData::SetTrialConstraintMultipliersFromPtr(
//     const SmartPtr<const Vector>& y_cptr,
//     const SmartPtr<const Vector>& y_dptr)
//   {
//     DBG_ASSERT(have_prototypes_);

//     trial_y_c_ = y_cptr;
//     trial_y_d_ = y_dptr;
//   }

//   void IpoptData::SetTrialBoundMultipliersFromPtr(
//     const SmartPtr<const Vector>& z_Lptr,
//     const SmartPtr<const Vector>& z_Uptr,
//     const SmartPtr<const Vector>& v_Lptr,
//     const SmartPtr<const Vector>& v_Uptr)
//   {
//     DBG_ASSERT(have_prototypes_);

//     trial_z_L_ = z_Lptr;
//     trial_z_U_ = z_Uptr;
//     trial_v_L_ = v_Lptr;
//     trial_v_U_ = v_Uptr;
//   }

  void IpoptData::CopyTrialToCurrent()
  {
    curr_ = trial_;
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
    delta_aff_x_ = NULL;
    delta_aff_s_ = NULL;
    delta_aff_y_c_ = NULL;
    delta_aff_y_d_ = NULL;
    delta_aff_z_L_ = NULL;
    delta_aff_z_U_ = NULL;
    delta_aff_v_L_ = NULL;
    delta_aff_v_U_ = NULL;

    have_deltas_ = false;
    have_affine_deltas_ = false;
  }

} // namespace Ipopt
