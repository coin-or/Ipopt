// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpIpoptCalculatedQuantities.hpp"
#ifdef OLD_C_HEADERS
#include <math.h>
#else
#include <cmath>
#endif

#include <limits>

// ToDo remove the following dependency
#include "IpRestoIpoptNLP.hpp"

namespace Ipopt
{
  DBG_SET_VERBOSITY(0);

  IpoptCalculatedQuantities::IpoptCalculatedQuantities
  (const SmartPtr<IpoptNLP>& ip_nlp,
   const SmartPtr<IpoptData>& ip_data)
      :
      ip_nlp_(ip_nlp),
      ip_data_(ip_data),

      curr_slack_x_L_cache_(1),
      curr_slack_x_U_cache_(1),
      curr_slack_s_L_cache_(1),
      curr_slack_s_U_cache_(1),
      trial_slack_x_L_cache_(1),
      trial_slack_x_U_cache_(1),
      trial_slack_s_L_cache_(1),
      trial_slack_s_U_cache_(1),
      num_adjusted_slack_x_L_(0),
      num_adjusted_slack_x_U_(0),
      num_adjusted_slack_s_L_(0),
      num_adjusted_slack_s_U_(0),

      curr_f_cache_(0),
      trial_f_cache_(0),
      curr_grad_f_cache_(2),

      curr_barrier_obj_cache_(2),
      trial_barrier_obj_cache_(5),
      curr_grad_barrier_obj_x_cache_(1),
      curr_grad_barrier_obj_s_cache_(1),

      curr_c_cache_(1),
      trial_c_cache_(2),
      curr_d_cache_(1),
      trial_d_cache_(2),
      curr_d_minus_s_cache_(1),
      trial_d_minus_s_cache_(1),
      curr_jac_c_cache_(1),
      curr_jac_d_cache_(1),
      curr_jac_cT_times_vec_cache_(2),
      curr_jac_dT_times_vec_cache_(2),
      curr_jac_c_times_vec_cache_(1),
      curr_jac_d_times_vec_cache_(1),
      curr_exact_hessian_cache_(1),
      curr_constraint_violation_cache_(2),
      trial_constraint_violation_cache_(5),

      curr_grad_lag_x_cache_(1),
      curr_grad_lag_s_cache_(1),
      curr_compl_x_L_cache_(1),
      curr_compl_x_U_cache_(1),
      curr_compl_s_L_cache_(1),
      curr_compl_s_U_cache_(1),
      curr_relaxed_compl_x_L_cache_(1),
      curr_relaxed_compl_x_U_cache_(1),
      curr_relaxed_compl_s_L_cache_(1),
      curr_relaxed_compl_s_U_cache_(1),
      curr_primal_infeasibility_cache_(3),
      trial_primal_infeasibility_cache_(3),
      curr_dual_infeasibility_cache_(3),
      curr_complementarity_cache_(6),
      curr_centrality_measure_cache_(1),
      curr_nlp_error_cache_(1),
      curr_barrier_error_cache_(1),
      curr_primal_dual_error_cache_(1),
      curr_relaxed_primal_dual_error_cache_(5),

      primal_frac_to_the_bound_cache_(5),
      dual_frac_to_the_bound_cache_(5),
      slack_frac_to_the_bound_cache_(5),

      curr_sigma_x_cache_(1),
      curr_sigma_s_cache_(1),

      curr_avrg_compl_cache_(1),
      trial_avrg_compl_cache_(1),
      curr_gradBarrTDelta_cache_(1),

      dampind_x_L_(NULL),
      dampind_x_U_(NULL),
      dampind_s_L_(NULL),
      dampind_s_U_(NULL),

      initialize_called_(false)
  {
    DBG_START_METH("IpoptCalculatedQuantities::IpoptCalculatedQuantities",
                   dbg_verbosity);
    DBG_ASSERT(IsValid(ip_nlp_) && IsValid(ip_data_));
  }

  IpoptCalculatedQuantities::~IpoptCalculatedQuantities()
  {}

  bool IpoptCalculatedQuantities::Initialize(const Journalist& jnlst,
      const OptionsList& options,
      const std::string& prefix)
  {
    Number value;
    std::string svalue;

    if (options.GetNumericValue("s_max", value, prefix)) {
      ASSERT_EXCEPTION(value > 0.0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"s_max\": This value must be larger than 0.");
      s_max_ = value;
    }
    else {
      s_max_ = 100.;
    }

    if (options.GetNumericValue("kappa_d", value, prefix)) {
      ASSERT_EXCEPTION(value >= 0.0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"kappa_d\": This value must be non-negative.");
      kappa_d_ = value;
    }
    else {
      kappa_d_ = 1e-5;
    }

    if (options.GetNumericValue("s_move_", value, prefix)) {
      ASSERT_EXCEPTION(value >= 0.0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"s_move\": This value must be non-negative.");
      s_move_ = value;
    }
    else {
      s_move_ = pow(std::numeric_limits<double>::epsilon(), 0.75);
    }

    if (options.GetValue("constraint_violation_normtype", svalue, prefix)) {
      if (svalue=="norm_1" || svalue=="1-norm" || svalue=="1norm") {
        constr_viol_normtype_ = NORM_1;
      }
      else if (svalue=="norm_2" || svalue=="2-norm" || svalue=="2norm") {
        constr_viol_normtype_ = NORM_2;
      }
      else if (svalue=="norm_max" || svalue=="max-norm" || svalue=="maxnorm") {
        constr_viol_normtype_ = NORM_MAX;
      }
      else {
        ASSERT_EXCEPTION(false, OptionsList::OPTION_OUT_OF_RANGE,
                         "Option \"constraint_violation_normtype\": Unknown keyword \""+svalue+"\".");
      }
    }
    else {
      constr_viol_normtype_ = NORM_1;
    }

    initialize_called_ = true;
    return true;
  }

  ///////////////////////////////////////////////////////////////////////////
  //                         Slack Calculations                            //
  ///////////////////////////////////////////////////////////////////////////

  SmartPtr<Vector>
  IpoptCalculatedQuantities::CalcSlack_L(const Matrix& P,
                                         const Vector& x,
                                         const Vector& x_bound)
  {
    DBG_START_METH("IpoptCalculatedQuantities::CalcSlack_L",
                   dbg_verbosity);
    SmartPtr<Vector> result;
    result = x_bound.MakeNew();
    result->Copy(x_bound);
    P.TransMultVector(1.0, x, -1.0, *result);
    return result;
  }

  SmartPtr<Vector>
  IpoptCalculatedQuantities::CalcSlack_U(const Matrix& P,
                                         const Vector& x,
                                         const Vector& x_bound)
  {
    DBG_START_METH("IpoptCalculatedQuantities::CalcSlack_U",
                   dbg_verbosity);
    SmartPtr<Vector> result;
    result = x_bound.MakeNew();
    result->Copy(x_bound);
    P.TransMultVector(-1.0, x, 1.0, *result);
    return result;
  }

  SmartPtr<const Vector> IpoptCalculatedQuantities::curr_slack_x_L()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_slack_x_L()",
                   dbg_verbosity);
    SmartPtr<Vector> result;
    SmartPtr<const Vector> x = ip_data_->curr_x();
    SmartPtr<const Vector> x_bound = ip_nlp_->x_L();
    if (!curr_slack_x_L_cache_.GetCachedResult1Dep(result, *x)) {
      if (!trial_slack_x_L_cache_.GetCachedResult1Dep(result, *x)) {
        SmartPtr<const Matrix> P = ip_nlp_->Px_L();
        DBG_PRINT_VECTOR(2,"x_L", *x_bound);
        result = CalcSlack_L(*P, *x, *x_bound);
        DBG_ASSERT(num_adjusted_slack_x_L_==0);
        num_adjusted_slack_x_L_ =
          CalculateSafeSlack(result, x_bound, x, ip_data_->curr_z_L());
      }
      curr_slack_x_L_cache_.AddCachedResult1Dep(result, *x);
    }
    return ConstPtr(result);
  }

  SmartPtr<const Vector> IpoptCalculatedQuantities::curr_slack_x_U()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_slack_x_U()",
                   dbg_verbosity);

    SmartPtr<Vector> result;
    SmartPtr<const Vector> x = ip_data_->curr_x();
    SmartPtr<const Vector> x_bound = ip_nlp_->x_U();
    if (!curr_slack_x_U_cache_.GetCachedResult1Dep(result, *x)) {
      if (!trial_slack_x_U_cache_.GetCachedResult1Dep(result, *x)) {
        SmartPtr<const Matrix> P = ip_nlp_->Px_U();
        result = CalcSlack_U(*P, *x, *x_bound);
        DBG_ASSERT(num_adjusted_slack_x_U_==0);
        num_adjusted_slack_x_U_ =
          CalculateSafeSlack(result, x_bound, x, ip_data_->curr_z_U());
      }
      curr_slack_x_U_cache_.AddCachedResult1Dep(result, *x);
    }
    return ConstPtr(result);
  }

  SmartPtr<const Vector> IpoptCalculatedQuantities::curr_slack_s_L()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_slack_s_L()",
                   dbg_verbosity);

    SmartPtr<Vector> result;
    SmartPtr<const Vector> s = ip_data_->curr_s();
    SmartPtr<const Vector> s_bound = ip_nlp_->d_L();
    if (!curr_slack_s_L_cache_.GetCachedResult1Dep(result, *s)) {
      if (!trial_slack_s_L_cache_.GetCachedResult1Dep(result, *s)) {
        SmartPtr<const Matrix> P = ip_nlp_->Pd_L();
        result = CalcSlack_L(*P, *s, *s_bound);
        DBG_ASSERT(num_adjusted_slack_s_L_==0);
        num_adjusted_slack_s_L_ =
          CalculateSafeSlack(result, s_bound, s, ip_data_->curr_v_L());
      }
      curr_slack_s_L_cache_.AddCachedResult1Dep(result, *s);
    }
    return ConstPtr(result);
  }

  SmartPtr<const Vector> IpoptCalculatedQuantities::curr_slack_s_U()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_slack_s_U()",
                   dbg_verbosity);

    SmartPtr<Vector> result;
    SmartPtr<const Vector> s = ip_data_->curr_s();
    SmartPtr<const Vector> s_bound = ip_nlp_->d_U();
    if (!curr_slack_s_U_cache_.GetCachedResult1Dep(result, *s)) {
      if (!trial_slack_s_U_cache_.GetCachedResult1Dep(result, *s)) {
        SmartPtr<const Matrix> P = ip_nlp_->Pd_U();
        DBG_PRINT_VECTOR(2, "s", *s);
        DBG_PRINT_VECTOR(2, "s_U", *s_bound);
        result = CalcSlack_U(*P, *s, *s_bound);
        DBG_ASSERT(num_adjusted_slack_s_U_==0);
        num_adjusted_slack_s_U_ =
          CalculateSafeSlack(result, s_bound, s, ip_data_->curr_v_U());
        DBG_PRINT_VECTOR(2, "result", *result);
      }
      curr_slack_s_U_cache_.AddCachedResult1Dep(result, *s);
    }
    return ConstPtr(result);
  }

  SmartPtr<const Vector> IpoptCalculatedQuantities::trial_slack_x_L()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_slack_x_L()",
                   dbg_verbosity);

    num_adjusted_slack_x_L_ = 0;
    SmartPtr<Vector> result;
    SmartPtr<const Vector> x = ip_data_->trial_x();
    SmartPtr<const Vector> x_bound = ip_nlp_->x_L();
    if (!trial_slack_x_L_cache_.GetCachedResult1Dep(result, *x)) {
      if (!curr_slack_x_L_cache_.GetCachedResult1Dep(result, *x)) {
        SmartPtr<const Matrix> P = ip_nlp_->Px_L();
        result = CalcSlack_L(*P, *x, *x_bound);
        DBG_ASSERT(num_adjusted_slack_x_L_==0);
        num_adjusted_slack_x_L_ =
          CalculateSafeSlack(result, x_bound, x, ip_data_->curr_z_L());
      }
      trial_slack_x_L_cache_.AddCachedResult1Dep(result, *x);
    }
    return ConstPtr(result);
  }

  SmartPtr<const Vector> IpoptCalculatedQuantities::trial_slack_x_U()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_slack_x_U()",
                   dbg_verbosity);

    num_adjusted_slack_x_U_ = 0;
    SmartPtr<Vector> result;
    SmartPtr<const Vector> x = ip_data_->trial_x();
    SmartPtr<const Vector> x_bound = ip_nlp_->x_U();
    if (!trial_slack_x_U_cache_.GetCachedResult1Dep(result, *x)) {
      if (!curr_slack_x_U_cache_.GetCachedResult1Dep(result, *x)) {
        SmartPtr<const Matrix> P = ip_nlp_->Px_U();
        result = CalcSlack_U(*P, *x, *x_bound);
        DBG_ASSERT(num_adjusted_slack_x_U_==0);
        num_adjusted_slack_x_U_ =
          CalculateSafeSlack(result, x_bound, x, ip_data_->curr_z_U());
      }
      trial_slack_x_U_cache_.AddCachedResult1Dep(result, *x);
    }
    return ConstPtr(result);
  }

  SmartPtr<const Vector> IpoptCalculatedQuantities::trial_slack_s_L()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_slack_s_L()",
                   dbg_verbosity);

    num_adjusted_slack_s_L_ = 0;
    SmartPtr<Vector> result;
    SmartPtr<const Vector> s = ip_data_->trial_s();
    SmartPtr<const Vector> s_bound = ip_nlp_->d_L();
    if (!trial_slack_s_L_cache_.GetCachedResult1Dep(result, *s)) {
      if (!curr_slack_s_L_cache_.GetCachedResult1Dep(result, *s)) {
        SmartPtr<const Matrix> P = ip_nlp_->Pd_L();
        result = CalcSlack_L(*P, *s, *s_bound);
        DBG_ASSERT(num_adjusted_slack_s_L_==0);
        num_adjusted_slack_s_L_ =
          CalculateSafeSlack(result, s_bound, s, ip_data_->curr_v_L());
      }
      trial_slack_s_L_cache_.AddCachedResult1Dep(result, *s);
    }
    return ConstPtr(result);
  }

  SmartPtr<const Vector> IpoptCalculatedQuantities::trial_slack_s_U()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_slack_s_U()",
                   dbg_verbosity);

    num_adjusted_slack_s_U_ = 0;
    SmartPtr<Vector> result;
    SmartPtr<const Vector> s = ip_data_->trial_s();
    SmartPtr<const Vector> s_bound = ip_nlp_->d_U();
    if (!trial_slack_s_U_cache_.GetCachedResult1Dep(result, *s)) {
      if (!curr_slack_s_U_cache_.GetCachedResult1Dep(result, *s)) {
        SmartPtr<const Matrix> P = ip_nlp_->Pd_U();
        result = CalcSlack_U(*P, *s, *s_bound);
        DBG_ASSERT(num_adjusted_slack_s_U_==0);
        num_adjusted_slack_s_U_ =
          CalculateSafeSlack(result, s_bound, s, ip_data_->curr_v_U());
      }
      trial_slack_s_U_cache_.AddCachedResult1Dep(result, *s);
    }
    return ConstPtr(result);
  }

  Index IpoptCalculatedQuantities::
  CalculateSafeSlack(SmartPtr<Vector>& slack,
                     const SmartPtr<const Vector>& bound,
                     const SmartPtr<const Vector>& curr_point,
                     const SmartPtr<const Vector>& multiplier)
  {
    DBG_START_METH("IpoptCalculatedQuantities::CalculateSafeSlack", dbg_verbosity);
    DBG_ASSERT(initialize_called_);
    Index retval = 0;
    if (slack->Dim() > 0) {
      Number min_slack = slack->Min();
      // TODO we need to make sure that this also works for non-monotone MUs
      Number s_min = std::numeric_limits<Number>::epsilon()
                     * Min(1., ip_data_->curr_mu());
      DBG_PRINT((1,"s_min = %g, min_slack=%g\n", s_min, min_slack));
      if (min_slack < s_min) {
        // Need to correct the slacks and calculate new bounds...
        SmartPtr<Vector> t = slack->MakeNew();
        t->Copy(*slack);
        t->AddScalar(-s_min);
        t->ElementWiseSgn();

        SmartPtr<Vector> zero_vec = t->MakeNew();
        zero_vec->Set(0.0);
        t->ElementWiseMin(*zero_vec);
        t->Scal(-1.0);
        retval = (Index)t->Asum();
        DBG_PRINT((1,"Number of slack corrections = %d\n", retval));
        DBG_PRINT_VECTOR(2, "t(sgn)", *t);

        SmartPtr<Vector> t2 = t->MakeNew();
        t2->Copy(*multiplier);
        t2->ElementWiseReciprocal();
        t2->Scal(ip_data_->curr_mu());

        SmartPtr<Vector> s_min_vec = t2->MakeNew();
        s_min_vec->Set(s_min);

        t2->ElementWiseMax(*s_min_vec);
        t2->Axpy(-1.0, *slack);
        DBG_PRINT_VECTOR(2, "tw(smin,mu/mult)", *t2);

        t->ElementWiseMultiply(*t2);
        t->Axpy(1.0, *slack);

        SmartPtr<Vector> t_max = t2;
        t_max->Set(1.0);
        t_max->ElementWiseMax(*bound);
        t_max->Scal(s_move_);
        t_max->Axpy(1.0, *slack);
        DBG_PRINT_VECTOR(2, "t_max", *t_max);

        t->ElementWiseMin(*t_max);
        DBG_PRINT_VECTOR(2, "new_slack", *t);

        slack = t;
        return retval;
      }
    }

    return retval;
  }

  Index
  IpoptCalculatedQuantities::AdjustedTrialSlacks()
  {
    DBG_START_METH("IpoptCalculatedQuantities::AdjustedTrialSlacks()",
		   dbg_verbosity);
    Index result =  (num_adjusted_slack_x_L_ +
		     num_adjusted_slack_x_U_ +
		     num_adjusted_slack_s_L_ +
		     num_adjusted_slack_s_U_);
    DBG_PRINT((1,"result = %d\n", result));
    return result;
  }

  void
  IpoptCalculatedQuantities::ResetAdjustedTrialSlacks()
  {
    num_adjusted_slack_x_L_
    = num_adjusted_slack_x_U_
      = num_adjusted_slack_s_L_
        = num_adjusted_slack_s_U_ = 0;
  }

  ///////////////////////////////////////////////////////////////////////////
  //                          Objective Function                           //
  ///////////////////////////////////////////////////////////////////////////

  Number
  IpoptCalculatedQuantities::curr_f()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_f()",
                   dbg_verbosity);
    Number result;
    SmartPtr<const Vector> x = ip_data_->curr_x();
    DBG_PRINT_VECTOR(2,"curr_x",*x);
    DBG_PRINT((1, "curr_x tag = %d\n", x->GetTag()));

    // ToDo: For now we make the value dependent on curr_mu during the
    // restoration phase (because Eta in the restoration phase
    // objective depends on it).  Need more elegant solution later.
    std::vector<const TaggedObject*> tdeps;
    tdeps.push_back(GetRawPtr(x));
    std::vector<Number> sdeps;
    if (in_restoration_phase()) {
      sdeps.push_back(ip_data_->curr_mu());
    }
    else {
      sdeps.push_back(-1.);
    }

    if (!curr_f_cache_.GetCachedResult(result, tdeps, sdeps)) {
      if (!trial_f_cache_.GetCachedResult(result, tdeps, sdeps)) {
	DBG_PRINT((2,"evaluate curr f\n"));
        result = ip_nlp_->f(*x);
      }
      curr_f_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    DBG_PRINT((1,"result (curr_f) = %e\n", result));
    return result;
  }

  Number
  IpoptCalculatedQuantities::trial_f()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_f()",
                   dbg_verbosity);
    Number result;
    SmartPtr<const Vector> x = ip_data_->trial_x();
    DBG_PRINT_VECTOR(2,"trial_x",*x);
    DBG_PRINT((1, "trial_x tag = %d\n", x->GetTag()));

    // ToDo: For now we make the value dependent on curr_mu during the
    // restoration phase (because Eta in the restoration phase
    // objective depends on it).  Need more elegant solution later.
    std::vector<const TaggedObject*> tdeps;
    tdeps.push_back(GetRawPtr(x));
    std::vector<Number> sdeps;
    if (in_restoration_phase()) {
      sdeps.push_back(ip_data_->curr_mu());
    }
    else {
      sdeps.push_back(-1.);
    }

    if (!trial_f_cache_.GetCachedResult(result, tdeps, sdeps)) {
      if (!curr_f_cache_.GetCachedResult(result, tdeps, sdeps)) {
	DBG_PRINT((2,"evaluate trial f\n"));
        result = ip_nlp_->f(*x);
      }
      trial_f_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    DBG_PRINT((1,"result (trial_f) = %e\n", result));
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_grad_f()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_grad_f()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->curr_x();

    // ToDo: For now we make the value dependent on curr_mu during the
    // restoration phase (because Eta in the restoration phase
    // objective depends on it).  Need more elegant solution later.
    std::vector<const TaggedObject*> tdeps;
    tdeps.push_back(GetRawPtr(x));
    std::vector<Number> sdeps;
    if (in_restoration_phase()) {
      sdeps.push_back(ip_data_->curr_mu());
    }
    else {
      sdeps.push_back(-1.);
    }

    if (!curr_grad_f_cache_.GetCachedResult(result, tdeps, sdeps)) {
      result = ip_nlp_->grad_f(*x);
      curr_grad_f_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    return result;
  }

  ///////////////////////////////////////////////////////////////////////////
  //                    Barrier Objective Function                         //
  ///////////////////////////////////////////////////////////////////////////
  Number
  IpoptCalculatedQuantities::CalcBarrierTerm(Number mu,
      const Vector& slack_x_L,
      const Vector& slack_x_U,
      const Vector& slack_s_L,
      const Vector& slack_s_U)
  {
    DBG_START_METH("IpoptCalculatedQuantities::CalcBarrierTerm",
                   dbg_verbosity);
    DBG_ASSERT(initialize_called_);

    DBG_PRINT_VECTOR(2, "slack_x_L", slack_x_L);
    DBG_PRINT_VECTOR(2, "slack_x_U", slack_x_U);
    DBG_PRINT_VECTOR(2, "slack_s_L", slack_s_L);
    DBG_PRINT_VECTOR(2, "slack_s_U", slack_s_U);
	
    Number retval=0.;
    retval += slack_x_L.SumLogs();
    retval += slack_x_U.SumLogs();
    retval += slack_s_L.SumLogs();
    retval += slack_s_U.SumLogs();
    retval *= -mu;

    // Include the linear damping term if kappa_d is nonzero.
    if (kappa_d_>0) {
      SmartPtr<const Vector> dampind_x_L;
      SmartPtr<const Vector> dampind_x_U;
      SmartPtr<const Vector> dampind_s_L;
      SmartPtr<const Vector> dampind_s_U;
      ComputeDampingIndicators(dampind_x_L, dampind_x_U, dampind_s_L, dampind_s_U);

      SmartPtr<Vector> tmp;
      tmp = slack_x_L.MakeNew();
      tmp->Copy(slack_x_L);
      tmp->ElementWiseMultiply(*dampind_x_L);
      retval += kappa_d_ * mu * tmp->Asum();
      tmp = slack_x_U.MakeNew();
      tmp->Copy(slack_x_U);
      tmp->ElementWiseMultiply(*dampind_x_U);
      retval += kappa_d_ * mu * tmp->Asum();
      tmp = slack_s_L.MakeNew();
      tmp->Copy(slack_s_L);
      tmp->ElementWiseMultiply(*dampind_s_L);
      retval += kappa_d_ * mu * tmp->Asum();
      tmp = slack_s_U.MakeNew();
      tmp->Copy(slack_s_U);
      tmp->ElementWiseMultiply(*dampind_s_U);
      retval += kappa_d_ * mu * tmp->Asum();
    }

    DBG_ASSERT(FiniteNumber(retval));
    return retval;
  }

  Number
  IpoptCalculatedQuantities::curr_barrier_obj()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_barrier_obj()",
                   dbg_verbosity);
    Number result;

    SmartPtr<const Vector> x = ip_data_->curr_x();
    SmartPtr<const Vector> s = ip_data_->curr_s();
    DBG_PRINT_VECTOR(2,"curr_x",*x);
    DBG_PRINT_VECTOR(2,"curr_s",*s);
    std::vector<const TaggedObject*> tdeps;
    tdeps.push_back(GetRawPtr(x));
    tdeps.push_back(GetRawPtr(s));

    Number mu = ip_data_->curr_mu();
    DBG_PRINT((1,"curr_mu=%e\n",mu));
    std::vector<Number> sdeps;
    sdeps.push_back(mu);

    if (!curr_barrier_obj_cache_.GetCachedResult(result, tdeps, sdeps)) {
      if (!trial_barrier_obj_cache_.GetCachedResult(result, tdeps, sdeps)) {
        result = curr_f();
	DBG_PRINT((1,"curr_F=%e\n",result));
        result += CalcBarrierTerm(mu,
                                  *curr_slack_x_L(),
                                  *curr_slack_x_U(),
                                  *curr_slack_s_L(),
                                  *curr_slack_s_U());
      }
      curr_barrier_obj_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    DBG_ASSERT(FiniteNumber(result));
    return result;
  }

  Number
  IpoptCalculatedQuantities::trial_barrier_obj()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_barrier_obj()",
                   dbg_verbosity);
    Number result;

    SmartPtr<const Vector> x = ip_data_->trial_x();
    SmartPtr<const Vector> s = ip_data_->trial_s();
    DBG_PRINT_VECTOR(2,"trial_x",*x);
    DBG_PRINT_VECTOR(2,"trial_s",*s);
    std::vector<const TaggedObject*> tdeps;
    tdeps.push_back(GetRawPtr(x));
    tdeps.push_back(GetRawPtr(s));

    Number mu = ip_data_->curr_mu();
    DBG_PRINT((1,"trial_mu=%e\n",mu));
    std::vector<Number> sdeps;
    sdeps.push_back(mu);

    if (!trial_barrier_obj_cache_.GetCachedResult(result, tdeps, sdeps)) {
      if (!curr_barrier_obj_cache_.GetCachedResult(result, tdeps, sdeps)) {
        result = trial_f();
	DBG_PRINT((1,"trial_F=%e\n",result));
        result += CalcBarrierTerm(ip_data_->curr_mu(),
                                  *trial_slack_x_L(),
                                  *trial_slack_x_U(),
                                  *trial_slack_s_L(),
                                  *trial_slack_s_U());
      }
      trial_barrier_obj_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    DBG_ASSERT(FiniteNumber(result));
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_grad_barrier_obj_x()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_grad_barrier_obj_x()",
                   dbg_verbosity);
    DBG_ASSERT(initialize_called_);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> x = ip_data_->curr_x();
    std::vector<const TaggedObject*> tdeps;
    tdeps.push_back(GetRawPtr(x));
    Number mu = ip_data_->curr_mu();
    std::vector<Number> sdeps;
    sdeps.push_back(mu);

    if (!curr_grad_barrier_obj_x_cache_.GetCachedResult(result, tdeps, sdeps)) {
      SmartPtr<Vector> tmp1 = x->MakeNew();
      tmp1->Copy(*curr_grad_f());

      SmartPtr<const Vector> slack = curr_slack_x_L();
      SmartPtr<Vector> tmp2 = slack->MakeNew();
      tmp2->Set(1.);
      tmp2->ElementWiseDivide(*slack);
      ip_nlp_->Px_L()->MultVector(-mu, *tmp2, 1., *tmp1);

      slack = curr_slack_x_U();
      tmp2 = slack->MakeNew();
      tmp2->Set(1.);
      tmp2->ElementWiseDivide(*slack);
      ip_nlp_->Px_U()->MultVector(mu, *tmp2, 1., *tmp1);

      // Take care of linear damping terms
      if (kappa_d_>0.) {
        SmartPtr<const Vector> dampind_x_L;
        SmartPtr<const Vector> dampind_x_U;
        SmartPtr<const Vector> dampind_s_L;
        SmartPtr<const Vector> dampind_s_U;
        ComputeDampingIndicators(dampind_x_L, dampind_x_U, dampind_s_L, dampind_s_U);

        ip_nlp_->Px_L()->MultVector(kappa_d_*mu, *dampind_x_L, 1., *tmp1);
        ip_nlp_->Px_U()->MultVector(-kappa_d_*mu, *dampind_x_U, 1., *tmp1);
      }

      result = ConstPtr(tmp1);

      curr_grad_barrier_obj_x_cache_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_grad_barrier_obj_s()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_grad_barrier_obj_s()",
                   dbg_verbosity);
    DBG_ASSERT(initialize_called_);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> s = ip_data_->curr_s();
    std::vector<const TaggedObject*> tdeps;
    tdeps.push_back(GetRawPtr(s));
    Number mu = ip_data_->curr_mu();
    std::vector<Number> sdeps;
    sdeps.push_back(mu);

    if (!curr_grad_barrier_obj_s_cache_.GetCachedResult(result, tdeps, sdeps)) {
      SmartPtr<Vector> tmp1 = s->MakeNew();

      SmartPtr<const Vector> slack = curr_slack_s_L();
      SmartPtr<Vector> tmp2 = slack->MakeNew();
      tmp2->Set(1.);
      tmp2->ElementWiseDivide(*slack);
      ip_nlp_->Pd_L()->MultVector(-mu, *tmp2, 0., *tmp1);

      slack = curr_slack_s_U();
      tmp2 = slack->MakeNew();
      tmp2->Set(1.);
      tmp2->ElementWiseDivide(*slack);
      ip_nlp_->Pd_U()->MultVector(mu, *tmp2, 1., *tmp1);

      // Take care of linear damping terms
      if (kappa_d_>0.) {
        SmartPtr<const Vector> dampind_x_L;
        SmartPtr<const Vector> dampind_x_U;
        SmartPtr<const Vector> dampind_s_L;
        SmartPtr<const Vector> dampind_s_U;
        ComputeDampingIndicators(dampind_x_L, dampind_x_U, dampind_s_L, dampind_s_U);

        ip_nlp_->Pd_L()->MultVector(kappa_d_*mu, *dampind_s_L, 1., *tmp1);
        ip_nlp_->Pd_U()->MultVector(-kappa_d_*mu, *dampind_s_U, 1., *tmp1);
      }

      result = ConstPtr(tmp1);

      curr_grad_barrier_obj_s_cache_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;
  }

  void
  IpoptCalculatedQuantities::ComputeDampingIndicators(SmartPtr<const Vector>& dampind_x_L,
      SmartPtr<const Vector>& dampind_x_U,
      SmartPtr<const Vector>& dampind_s_L,
      SmartPtr<const Vector>& dampind_s_U)
  {
    DBG_START_METH("IpoptCalculatedQuantities::ComputeDampingFilters()",
                   dbg_verbosity);

    // Assume that all indicators have to be computed if one of the
    // SmartPtrs is still zero.
    if (IsNull(dampind_x_L_)) {
      // First for x
      SmartPtr<Vector> diff = ip_data_->curr_x()->MakeNew();

      SmartPtr<Vector> vones = ip_nlp_->x_L()->MakeNew();
      vones->Set(1.0);
      ip_nlp_->Px_L()->MultVector(1.0, *vones, 0.0, *diff);
      vones = ip_nlp_->x_U()->MakeNew();
      vones->Set(1.0);
      ip_nlp_->Px_U()->MultVector(-1.0, *vones, 1.0, *diff);

      dampind_x_L_ = ip_nlp_->x_L()->MakeNew();
      ip_nlp_->Px_L()->TransMultVector(1.0, *diff, 0.0, *dampind_x_L_);

      dampind_x_U_ = ip_nlp_->x_U()->MakeNew();
      ip_nlp_->Px_U()->TransMultVector(-1.0, *diff, 0.0, *dampind_x_U_);

      // New for s
      diff = ip_data_->curr_s()->MakeNew();

      vones = ip_nlp_->d_L()->MakeNew();
      vones->Set(1.0);
      ip_nlp_->Pd_L()->MultVector(1.0, *vones, 0.0, *diff);
      vones = ip_nlp_->d_U()->MakeNew();
      vones->Set(1.0);
      ip_nlp_->Pd_U()->MultVector(-1.0, *vones, 1.0, *diff);

      dampind_s_L_ = ip_nlp_->d_L()->MakeNew();
      ip_nlp_->Pd_L()->TransMultVector(1.0, *diff, 0.0, *dampind_s_L_);

      dampind_s_U_ = ip_nlp_->d_U()->MakeNew();
      ip_nlp_->Pd_U()->TransMultVector(-1.0, *diff, 0.0, *dampind_s_U_);

      DBG_PRINT_VECTOR(2, "dampind_x_L_", *dampind_x_L_);
      DBG_PRINT_VECTOR(2, "dampind_x_U_", *dampind_x_U_);
      DBG_PRINT_VECTOR(2, "dampind_s_L_", *dampind_s_L_);
      DBG_PRINT_VECTOR(2, "dampind_s_U_", *dampind_s_U_);
    }

    dampind_x_L = ConstPtr(dampind_x_L_);
    dampind_x_U = ConstPtr(dampind_x_U_);
    dampind_s_L = ConstPtr(dampind_s_L_);
    dampind_s_U = ConstPtr(dampind_s_U_);
  }

  ///////////////////////////////////////////////////////////////////////////
  //                                Constraints                            //
  ///////////////////////////////////////////////////////////////////////////

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_c()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_c()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->curr_x();

    if (!curr_c_cache_.GetCachedResult1Dep(result, *x)) {
      if (!trial_c_cache_.GetCachedResult1Dep(result, *x)) {
        result = ip_nlp_->c(*x);
      }
      curr_c_cache_.AddCachedResult1Dep(result, *x);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::trial_c()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_c()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->trial_x();

    if (!trial_c_cache_.GetCachedResult1Dep(result, *x)) {
      if (!curr_c_cache_.GetCachedResult1Dep(result, *x)) {
        result = ip_nlp_->c(*x);
      }
      trial_c_cache_.AddCachedResult1Dep(result, *x);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_d()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_d()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->curr_x();

    if (!curr_d_cache_.GetCachedResult1Dep(result, *x)) {
      if (!trial_d_cache_.GetCachedResult1Dep(result, *x)) {
        result = ip_nlp_->d(*x);
      }
      curr_d_cache_.AddCachedResult1Dep(result, *x);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::trial_d()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_d()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->trial_x();

    if (!trial_d_cache_.GetCachedResult1Dep(result, *x)) {
      if (!curr_d_cache_.GetCachedResult1Dep(result, *x)) {
        result = ip_nlp_->d(*x);
      }
      trial_d_cache_.AddCachedResult1Dep(result, *x);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_d_minus_s()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_d_minus_s()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> x = ip_data_->curr_x();
    SmartPtr<const Vector> s = ip_data_->curr_s();

    if (!curr_d_minus_s_cache_.GetCachedResult2Dep(result, *x, *s)) {
      if (!trial_d_minus_s_cache_.GetCachedResult2Dep(result, *x, *s)) {
        SmartPtr<Vector> tmp = s->MakeNew();
        tmp->Copy(*curr_d());
        tmp->Axpy(-1., *s);
        result = ConstPtr(tmp);
      }
      curr_d_minus_s_cache_.AddCachedResult2Dep(result, *x, *s);
    }

    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::trial_d_minus_s()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_d_minus_s()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> x = ip_data_->trial_x();
    SmartPtr<const Vector> s = ip_data_->trial_s();

    if (!trial_d_minus_s_cache_.GetCachedResult2Dep(result, *x, *s)) {
      if (!curr_d_minus_s_cache_.GetCachedResult2Dep(result, *x, *s)) {
        SmartPtr<Vector> tmp = s->MakeNew();
        tmp->Copy(*trial_d());
        tmp->Axpy(-1., *s);
        result = ConstPtr(tmp);
      }
      trial_d_minus_s_cache_.AddCachedResult2Dep(result, *x, *s);
    }

    return result;
  }

  SmartPtr<const Matrix>
  IpoptCalculatedQuantities::curr_jac_c()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_jac_c()",
                   dbg_verbosity);
    SmartPtr<const Matrix> result;
    SmartPtr<const Vector> x = ip_data_->curr_x();

    if (!curr_jac_c_cache_.GetCachedResult1Dep(result, *x)) {
      result = ip_nlp_->jac_c(*x);
      curr_jac_c_cache_.AddCachedResult1Dep(result, *x);
    }
    return result;
  }

  SmartPtr<const Matrix>
  IpoptCalculatedQuantities::curr_jac_d()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_jac_d()",
                   dbg_verbosity);
    SmartPtr<const Matrix> result;
    SmartPtr<const Vector> x = ip_data_->curr_x();

    if (!curr_jac_d_cache_.GetCachedResult1Dep(result, *x)) {
      result = ip_nlp_->jac_d(*x);
      curr_jac_d_cache_.AddCachedResult1Dep(result, *x);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_jac_c_times_vec(const Vector& vec)
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_jac_c_times_vec",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->curr_x();

    if (!curr_jac_c_times_vec_cache_.GetCachedResult2Dep(result, *x, vec)) {
      SmartPtr<Vector> tmp = ip_data_->curr_y_c()->MakeNew();
      curr_jac_c()->MultVector(1.0, vec, 0., *tmp);
      result = ConstPtr(tmp);
      curr_jac_c_times_vec_cache_.AddCachedResult2Dep(result, *x, vec);
    }

    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_jac_d_times_vec(const Vector& vec)
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_jac_d_times_vec()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->curr_x();

    if (!curr_jac_d_times_vec_cache_.GetCachedResult2Dep(result, *x, vec)) {
      SmartPtr<Vector> tmp = ip_data_->curr_s()->MakeNew();
      curr_jac_d()->MultVector(1.0, vec, 0., *tmp);
      result = ConstPtr(tmp);
      curr_jac_d_times_vec_cache_.AddCachedResult2Dep(result, *x, vec);
    }

    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_jac_cT_times_curr_y_c()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_jac_cT_times_curr_y_c()",
                   dbg_verbosity);
    return curr_jac_cT_times_vec(*ip_data_->curr_y_c());
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_jac_dT_times_curr_y_d()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_jac_dT_times_curr_y_d()",
                   dbg_verbosity);
    return curr_jac_dT_times_vec(*ip_data_->curr_y_d());
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_jac_cT_times_vec(const Vector& vec)
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_jac_cT_times_vec",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->curr_x();

    if (!curr_jac_cT_times_vec_cache_.GetCachedResult2Dep(result, *x, vec)) {
      SmartPtr<Vector> tmp = ip_data_->curr_x()->MakeNew();
      curr_jac_c()->TransMultVector(1.0, vec, 0., *tmp);
      result = ConstPtr(tmp);
      curr_jac_cT_times_vec_cache_.AddCachedResult2Dep(result, *x, vec);
    }

    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_jac_dT_times_vec(const Vector& vec)
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_jac_dT_times_vec()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->curr_x();

    if (!curr_jac_dT_times_vec_cache_.GetCachedResult2Dep(result, *x, vec)) {
      SmartPtr<Vector> tmp = ip_data_->curr_x()->MakeNew();
      curr_jac_d()->TransMultVector(1.0, vec, 0., *tmp);
      result = ConstPtr(tmp);
      curr_jac_dT_times_vec_cache_.AddCachedResult2Dep(result, *x, vec);
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::curr_constraint_violation()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_constraint_violation()",
                   dbg_verbosity);
    return curr_primal_infeasibility(constr_viol_normtype_);
  }

  Number
  IpoptCalculatedQuantities::trial_constraint_violation()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_constraint_violation()",
                   dbg_verbosity);
    return trial_primal_infeasibility(constr_viol_normtype_);
  }

  ///////////////////////////////////////////////////////////////////////////
  //                Exact Hessian using second derivatives                 //
  ///////////////////////////////////////////////////////////////////////////

  SmartPtr<const SymMatrix>
  IpoptCalculatedQuantities::curr_exact_hessian()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_exact_hessian()",
                   dbg_verbosity);

    SmartPtr<const SymMatrix> result;
    SmartPtr<const Vector> x = ip_data_->curr_x();
    SmartPtr<const Vector> y_c = ip_data_->curr_y_c();
    SmartPtr<const Vector> y_d = ip_data_->curr_y_d();

    // ToDo: For now we make the value dependent on curr_mu during the
    // restoration phase (because Eta in the restoration phase
    // objective depends on it).  Need more elegant solution later.
    std::vector<const TaggedObject*> tdeps;
    tdeps.push_back(GetRawPtr(x));
    tdeps.push_back(GetRawPtr(y_c));
    tdeps.push_back(GetRawPtr(y_d));
    std::vector<Number> sdeps;
    if (in_restoration_phase()) {
      sdeps.push_back(ip_data_->curr_mu());
    }
    else {
      sdeps.push_back(-1.);
    }

    if (!curr_exact_hessian_cache_.GetCachedResult(result, tdeps, sdeps)) {
      result = ip_nlp_->h(*x, 1.0, *y_c, *y_d);
      curr_exact_hessian_cache_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;
  }

  ///////////////////////////////////////////////////////////////////////////
  //                             Zero Hessian                              //
  ///////////////////////////////////////////////////////////////////////////

  SmartPtr<const SymMatrix>
  IpoptCalculatedQuantities::zero_hessian()
  {
    DBG_START_METH("IpoptCalculatedQuantities::zero_hessian()",
                   dbg_verbosity);

    SmartPtr<const Vector> x = ip_data_->curr_x();
    SmartPtr<Vector> zero_y_c = ip_data_->curr_y_c()->MakeNew();
    SmartPtr<Vector> zero_y_d = ip_data_->curr_y_d()->MakeNew();
    zero_y_c->Set(0.);
    zero_y_d->Set(0.);

    SmartPtr<const SymMatrix> h = ip_nlp_->h(*x, 0.0, *zero_y_c, *zero_y_d);

    DBG_PRINT_MATRIX(2, "zero_hessian", *h);

    return h;
  }

  ///////////////////////////////////////////////////////////////////////////
  //                  Optimality Error and its components                  //
  ///////////////////////////////////////////////////////////////////////////

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_grad_lag_x()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_grad_lag_x()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> x = ip_data_->curr_x();
    SmartPtr<const Vector> y_c = ip_data_->curr_y_c();
    SmartPtr<const Vector> y_d = ip_data_->curr_y_d();
    SmartPtr<const Vector> z_L = ip_data_->curr_z_L();
    SmartPtr<const Vector> z_U = ip_data_->curr_z_U();

    std::vector<const TaggedObject*> deps;
    deps.push_back(GetRawPtr(x));
    deps.push_back(GetRawPtr(y_c));
    deps.push_back(GetRawPtr(y_d));
    deps.push_back(GetRawPtr(z_L));
    deps.push_back(GetRawPtr(z_U));

    if (!curr_grad_lag_x_cache_.GetCachedResult(result, deps)) {
      SmartPtr<Vector> tmp = x->MakeNew();
      DBG_PRINT_VECTOR(2,"curr_grad_f",*curr_grad_f());
      tmp->Copy(*curr_grad_f());
      tmp->Axpy(1., *curr_jac_cT_times_curr_y_c());
      tmp->Axpy(1., *curr_jac_dT_times_curr_y_d());
      DBG_PRINT_VECTOR(2,"jac_cT*y_c",*curr_jac_cT_times_curr_y_c());
      DBG_PRINT_VECTOR(2,"jac_dT*y_d",*curr_jac_dT_times_curr_y_d());
      ip_nlp_->Px_L()->MultVector(-1., *ip_data_->curr_z_L(), 1., *tmp);
      ip_nlp_->Px_U()->MultVector(1., *ip_data_->curr_z_U(), 1., *tmp);
      result = ConstPtr(tmp);
      curr_grad_lag_x_cache_.AddCachedResult(result, deps);
    }

    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_grad_lag_s()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_grad_lag_s()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> y_d = ip_data_->curr_y_d();
    SmartPtr<const Vector> v_L = ip_data_->curr_v_L();
    SmartPtr<const Vector> v_U = ip_data_->curr_v_U();

    std::vector<const TaggedObject*> deps;
    deps.push_back(GetRawPtr(y_d));
    deps.push_back(GetRawPtr(v_L));
    deps.push_back(GetRawPtr(v_U));

    if (!curr_grad_lag_s_cache_.GetCachedResult(result, deps)) {
      SmartPtr<Vector> tmp = y_d->MakeNew();
      ip_nlp_->Pd_U()->MultVector(1., *ip_data_->curr_v_U(), 0., *tmp);
      ip_nlp_->Pd_L()->MultVector(-1., *ip_data_->curr_v_L(), 1., *tmp);
      tmp->Axpy(-1., *y_d);
      result = ConstPtr(tmp);
      curr_grad_lag_s_cache_.AddCachedResult(result, deps);
    }

    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::CalcCompl(const Vector& slack,
                                       const Vector& mult)
  {
    DBG_START_METH("IpoptCalculatedQuantities::CalcCompl()",
                   dbg_verbosity);
    SmartPtr<Vector> result = slack.MakeNew();
    result->Copy(slack);
    result->ElementWiseMultiply(mult);
    return ConstPtr(result);
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_compl_x_L()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_compl_x_L()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> slack = curr_slack_x_L();
    SmartPtr<const Vector> mult = ip_data_->curr_z_L();
    DBG_PRINT_VECTOR(2, "slack_x_L", *slack);
    DBG_PRINT_VECTOR(2, "z_L", *mult);

    if (!curr_compl_x_L_cache_.GetCachedResult2Dep(result,
        *slack, *mult)) {
      result = CalcCompl(*slack, *mult);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_compl_x_U()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_compl_x_U()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> slack = curr_slack_x_U();
    SmartPtr<const Vector> mult = ip_data_->curr_z_U();

    if (!curr_compl_x_U_cache_.GetCachedResult2Dep(result,
        *slack, *mult)) {
      result = CalcCompl(*slack, *mult);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_compl_s_L()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_compl_s_L()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> slack = curr_slack_s_L();
    SmartPtr<const Vector> mult = ip_data_->curr_v_L();

    if (!curr_compl_s_L_cache_.GetCachedResult2Dep(result,
        *slack, *mult)) {
      result = CalcCompl(*slack, *mult);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_compl_s_U()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_compl_s_U()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> slack = curr_slack_s_U();
    SmartPtr<const Vector> mult = ip_data_->curr_v_U();

    if (!curr_compl_s_U_cache_.GetCachedResult2Dep(result,
        *slack, *mult)) {
      DBG_PRINT_VECTOR(2,"slack_s_U", *slack);
      DBG_PRINT_VECTOR(2,"v_U", *mult);
      result = CalcCompl(*slack, *mult);
      DBG_PRINT_VECTOR(2,"result", *result);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_relaxed_compl_x_L()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_relaxed_compl_x_L()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> slack = curr_slack_x_L();
    SmartPtr<const Vector> mult = ip_data_->curr_z_L();
    std::vector<const TaggedObject*> tdeps;
    tdeps.push_back(GetRawPtr(slack));
    tdeps.push_back(GetRawPtr(mult));

    Number mu = ip_data_->curr_mu();
    std::vector<Number> sdeps;
    sdeps.push_back(mu);

    if (!curr_relaxed_compl_x_L_cache_.GetCachedResult(result, tdeps, sdeps)) {
      SmartPtr<Vector> tmp = slack->MakeNew();
      tmp->Copy(*curr_compl_x_L());
      tmp->AddScalar(-mu);
      result = ConstPtr(tmp);
      curr_relaxed_compl_x_L_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_relaxed_compl_x_U()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_relaxed_compl_x_U()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> slack = curr_slack_x_U();
    SmartPtr<const Vector> mult = ip_data_->curr_z_U();
    std::vector<const TaggedObject*> tdeps;
    tdeps.push_back(GetRawPtr(slack));
    tdeps.push_back(GetRawPtr(mult));

    Number mu = ip_data_->curr_mu();
    std::vector<Number> sdeps;
    sdeps.push_back(mu);

    if (!curr_relaxed_compl_x_U_cache_.GetCachedResult(result, tdeps, sdeps)) {
      SmartPtr<Vector> tmp = slack->MakeNew();
      tmp->Copy(*curr_compl_x_U());
      tmp->AddScalar(-mu);
      result = ConstPtr(tmp);
      curr_relaxed_compl_x_U_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_relaxed_compl_s_L()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_relaxed_compl_s_L()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> slack = curr_slack_s_L();
    SmartPtr<const Vector> mult = ip_data_->curr_v_L();
    std::vector<const TaggedObject*> tdeps;
    tdeps.push_back(GetRawPtr(slack));
    tdeps.push_back(GetRawPtr(mult));

    Number mu = ip_data_->curr_mu();
    std::vector<Number> sdeps;
    sdeps.push_back(mu);

    if (!curr_relaxed_compl_s_L_cache_.GetCachedResult(result, tdeps, sdeps)) {
      SmartPtr<Vector> tmp = slack->MakeNew();
      tmp->Copy(*curr_compl_s_L());
      tmp->AddScalar(-mu);
      result = ConstPtr(tmp);
      curr_relaxed_compl_s_L_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_relaxed_compl_s_U()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_relaxed_compl_s_U()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> slack = curr_slack_s_U();
    SmartPtr<const Vector> mult = ip_data_->curr_v_U();
    std::vector<const TaggedObject*> tdeps;
    tdeps.push_back(GetRawPtr(slack));
    tdeps.push_back(GetRawPtr(mult));

    Number mu = ip_data_->curr_mu();
    std::vector<Number> sdeps;
    sdeps.push_back(mu);

    if (!curr_relaxed_compl_s_U_cache_.GetCachedResult(result, tdeps, sdeps)) {
      SmartPtr<Vector> tmp = slack->MakeNew();
      tmp->Copy(*curr_compl_s_U());
      tmp->AddScalar(-mu);
      result = ConstPtr(tmp);
      curr_relaxed_compl_s_U_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    return result;
  }

  Number
  IpoptCalculatedQuantities::CalcNormOfType
  (IpoptCalculatedQuantities::ENormType NormType,
   const Vector& vec1, const Vector& vec2)
  {
    std::vector<SmartPtr<const Vector> > vecs;
    vecs.push_back(&vec1);
    vecs.push_back(&vec2);

    return CalcNormOfType(NormType, vecs);
  }

  Number
  IpoptCalculatedQuantities::CalcNormOfType
  (IpoptCalculatedQuantities::ENormType NormType,
   std::vector<SmartPtr<const Vector> > vecs)
  {
    Number result=0.;

    switch (NormType) {
      case NORM_1 :
      for (Index i=0; i<(Index)vecs.size(); i++) {
        result += vecs[i]->Asum();
      }
      break;
      case NORM_2 :
      for (Index i=0; i<(Index)vecs.size(); i++) {
        Number nrm = vecs[i]->Nrm2();
        result += nrm*nrm;
      }
      result = sqrt(result);
      break;
      case NORM_MAX :
      for (Index i=0; i<(Index)vecs.size(); i++) {
        result = Max(result, vecs[i]->Amax());
      }
      break;
      default:
      DBG_ASSERT("Unknown NormType.");
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::curr_primal_infeasibility
  (IpoptCalculatedQuantities::ENormType NormType)
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_primal_infeasibility()",
                   dbg_verbosity);
    Number result;

    SmartPtr<const Vector> x = ip_data_->curr_x();
    SmartPtr<const Vector> s = ip_data_->curr_s();

    DBG_PRINT_VECTOR(2, "x to eval", *x);
    DBG_PRINT_VECTOR(2, "s to eval", *s);
    DBG_PRINT((1,"NormType = %d\n", NormType))

    std::vector<const TaggedObject*> deps;
    deps.push_back(GetRawPtr(x));
    deps.push_back(GetRawPtr(s));
    std::vector<Number> sdeps;
    sdeps.push_back((Number)NormType);

    if (!curr_primal_infeasibility_cache_.GetCachedResult(result, deps, sdeps)) {
      if (!trial_primal_infeasibility_cache_.GetCachedResult(result, deps, sdeps)) {
        DBG_PRINT((1,"Recomputing recomputing infeasibility.\n"));
        SmartPtr<const Vector> c = curr_c();
        SmartPtr<const Vector> d_minus_s = curr_d_minus_s();

        DBG_PRINT_VECTOR(2,"c", *c);
        DBG_PRINT_VECTOR(2,"d_minus_s", *d_minus_s);

        result = CalcNormOfType(NormType, *c, *d_minus_s);

      }
      curr_primal_infeasibility_cache_.AddCachedResult(result, deps, sdeps);
    }

    DBG_PRINT((1,"result = %e\n",result));
    return result;
  }

  Number
  IpoptCalculatedQuantities::trial_primal_infeasibility
  (IpoptCalculatedQuantities::ENormType NormType)
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_primal_infeasibility()",
                   dbg_verbosity);
    Number result;

    SmartPtr<const Vector> x = ip_data_->trial_x();
    SmartPtr<const Vector> s = ip_data_->trial_s();

    DBG_PRINT_VECTOR(2, "x to eval", *x);
    DBG_PRINT_VECTOR(2, "s to eval", *s);
    DBG_PRINT((1,"NormType = %d\n", NormType))

    std::vector<const TaggedObject*> deps;
    deps.push_back(GetRawPtr(x));
    deps.push_back(GetRawPtr(s));
    std::vector<Number> sdeps;
    sdeps.push_back((Number)NormType);

    if (!trial_primal_infeasibility_cache_.GetCachedResult(result, deps, sdeps)) {
      if (!curr_primal_infeasibility_cache_.GetCachedResult(result, deps, sdeps)) {
        DBG_PRINT((1,"Recomputing recomputing infeasibility.\n"));
        SmartPtr<const Vector> c = trial_c();
        SmartPtr<const Vector> d_minus_s = trial_d_minus_s();

        DBG_PRINT_VECTOR(2,"c", *c);
        DBG_PRINT_VECTOR(2,"d_minus_s", *d_minus_s);

        result = CalcNormOfType(NormType, *c, *d_minus_s);
      }
      trial_primal_infeasibility_cache_.AddCachedResult(result, deps, sdeps);
    }

    DBG_PRINT((1,"result = %e\n",result));
    return result;
  }

  Number
  IpoptCalculatedQuantities::curr_dual_infeasibility
  (IpoptCalculatedQuantities::ENormType NormType)
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_dual_infeasibility()",
                   dbg_verbosity);
    Number result;

    SmartPtr<const Vector> x = ip_data_->curr_x();
    SmartPtr<const Vector> s = ip_data_->curr_s();
    SmartPtr<const Vector> y_c = ip_data_->curr_y_c();
    SmartPtr<const Vector> y_d = ip_data_->curr_y_d();
    SmartPtr<const Vector> z_L = ip_data_->curr_z_L();
    SmartPtr<const Vector> z_U = ip_data_->curr_z_U();
    SmartPtr<const Vector> v_L = ip_data_->curr_v_L();
    SmartPtr<const Vector> v_U = ip_data_->curr_v_U();

    std::vector<const TaggedObject*> deps;
    deps.push_back(GetRawPtr(x));
    deps.push_back(GetRawPtr(s));
    deps.push_back(GetRawPtr(y_c));
    deps.push_back(GetRawPtr(y_c));
    deps.push_back(GetRawPtr(z_L));
    deps.push_back(GetRawPtr(z_U));
    deps.push_back(GetRawPtr(v_L));
    deps.push_back(GetRawPtr(v_U));
    std::vector<Number> sdeps;
    sdeps.push_back((Number)NormType);

    if (!curr_dual_infeasibility_cache_.GetCachedResult(result, deps, sdeps)) {
      SmartPtr<const Vector> grad_lag_x = curr_grad_lag_x();
      SmartPtr<const Vector> grad_lag_s = curr_grad_lag_s();

      result = CalcNormOfType(NormType, *grad_lag_x, *grad_lag_s);

      curr_dual_infeasibility_cache_.AddCachedResult(result, deps, sdeps);
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::curr_complementarity
  (Number mu, IpoptCalculatedQuantities::ENormType NormType)
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_complementarity()",
                   dbg_verbosity);
    Number result;

    SmartPtr<const Vector> x = ip_data_->curr_x();
    SmartPtr<const Vector> s = ip_data_->curr_s();
    SmartPtr<const Vector> z_L = ip_data_->curr_z_L();
    SmartPtr<const Vector> z_U = ip_data_->curr_z_U();
    SmartPtr<const Vector> v_L = ip_data_->curr_v_L();
    SmartPtr<const Vector> v_U = ip_data_->curr_v_U();

    std::vector<const TaggedObject*> deps;
    deps.push_back(GetRawPtr(x));
    deps.push_back(GetRawPtr(s));
    deps.push_back(GetRawPtr(z_L));
    deps.push_back(GetRawPtr(z_U));
    deps.push_back(GetRawPtr(v_L));
    deps.push_back(GetRawPtr(v_U));
    std::vector<Number> sdeps;
    sdeps.push_back((Number)NormType);
    sdeps.push_back(mu);

    if (!curr_complementarity_cache_.GetCachedResult(result, deps, sdeps)) {

      std::vector<SmartPtr<const Vector> > vecs;
      SmartPtr<const Vector> compl_x_L = curr_compl_x_L();
      SmartPtr<const Vector> compl_x_U = curr_compl_x_U();
      SmartPtr<const Vector> compl_s_L = curr_compl_s_L();
      SmartPtr<const Vector> compl_s_U = curr_compl_s_U();

      if (mu==.0) {
        vecs.push_back(GetRawPtr(compl_x_L));
        vecs.push_back(GetRawPtr(compl_x_U));
        vecs.push_back(GetRawPtr(compl_s_L));
        vecs.push_back(GetRawPtr(compl_s_U));
      }
      else {
        SmartPtr<Vector> tmp = compl_x_L->MakeNew();
        tmp->Copy(*compl_x_L);
        tmp->AddScalar(-mu);
        vecs.push_back(GetRawPtr(tmp));
        tmp = compl_x_U->MakeNew();
        tmp->Copy(*compl_x_U);
        tmp->AddScalar(-mu);
        vecs.push_back(GetRawPtr(tmp));
        tmp = compl_s_L->MakeNew();
        tmp->Copy(*compl_s_L);
        tmp->AddScalar(-mu);
        vecs.push_back(GetRawPtr(tmp));
        tmp = compl_s_U->MakeNew();
        tmp->Copy(*compl_s_U);
        tmp->AddScalar(-mu);
        vecs.push_back(GetRawPtr(tmp));
      }

      result = CalcNormOfType(NormType, vecs);

      curr_complementarity_cache_.AddCachedResult(result, deps, sdeps);
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::CalcCentralityMeasure(const Vector& compl_x_L,
      const Vector& compl_x_U,
      const Vector& compl_s_L,
      const Vector& compl_s_U)
  {
    DBG_START_METH("IpoptCalculatedQuantities::CalcCentralityMeasure()",
                   dbg_verbosity);

    Number MinCompl = std::numeric_limits<Number>::max();
    bool have_bounds = false;

    Index n_compl_x_L = compl_x_L.Dim();
    Index n_compl_x_U = compl_x_U.Dim();
    Index n_compl_s_L = compl_s_L.Dim();
    Index n_compl_s_U = compl_s_U.Dim();

    // Compute the Minimum of all complementarities
    if( n_compl_x_L>0 ) {
      if( have_bounds ) {
        MinCompl = Min(MinCompl, compl_x_L.Min());
      }
      else {
        MinCompl = compl_x_L.Min();
      }
      have_bounds = true;
    }
    if( n_compl_x_U>0 ) {
      if( have_bounds ) {
        MinCompl = Min(MinCompl, compl_x_U.Min());
      }
      else {
        MinCompl = compl_x_U.Min();
      }
      have_bounds = true;
    }
    if( n_compl_s_L>0 ) {
      if( have_bounds ) {
        MinCompl = Min(MinCompl, compl_s_L.Min());
      }
      else {
        MinCompl = compl_s_L.Min();
      }
      have_bounds = true;
    }
    if( n_compl_s_U>0 ) {
      if( have_bounds ) {
        MinCompl = Min(MinCompl, compl_s_U.Min());
      }
      else {
        MinCompl = compl_s_U.Min();
      }
      have_bounds = true;
    }

    // If there are no bounds, just return 0.;
    if (!have_bounds) {
      return 0.;
    }

    DBG_PRINT_VECTOR(2, "compl_x_L", compl_x_L);
    DBG_PRINT_VECTOR(2, "compl_x_U", compl_x_U);
    DBG_PRINT_VECTOR(2, "compl_s_L", compl_s_L);
    DBG_PRINT_VECTOR(2, "compl_s_U", compl_s_U);

    DBG_ASSERT(MinCompl>0. && "There is a zero complementarity entry");

    Number avrg_compl = (compl_x_L.Asum() + compl_x_U.Asum() +
                         compl_s_L.Asum() + compl_s_U.Asum());
    DBG_PRINT((1,"sum_compl = %25.16e\n", avrg_compl));
    avrg_compl /= (n_compl_x_L + n_compl_x_U + n_compl_s_L + n_compl_s_U);
    DBG_PRINT((1,"avrg_compl = %25.16e\n", avrg_compl));
    DBG_PRINT((1,"MinCompl = %25.16e\n", MinCompl));

    Number xi = MinCompl/avrg_compl;
    // The folloking line added for the case that avrg_compl is
    // slighly smaller than MinCompl, due to numerical roundoff
    xi = Min(1., xi);

    return xi;
  }

  Number
  IpoptCalculatedQuantities::curr_centrality_measure()
  {
    Number result;

    SmartPtr<const Vector> x = ip_data_->curr_x();
    SmartPtr<const Vector> s = ip_data_->curr_s();
    SmartPtr<const Vector> z_L = ip_data_->curr_z_L();
    SmartPtr<const Vector> z_U = ip_data_->curr_z_U();
    SmartPtr<const Vector> v_L = ip_data_->curr_v_L();
    SmartPtr<const Vector> v_U = ip_data_->curr_v_U();

    std::vector<const TaggedObject*> tdeps;
    tdeps.push_back(GetRawPtr(x));
    tdeps.push_back(GetRawPtr(s));
    tdeps.push_back(GetRawPtr(z_L));
    tdeps.push_back(GetRawPtr(z_U));
    tdeps.push_back(GetRawPtr(v_L));
    tdeps.push_back(GetRawPtr(z_U));

    if (!curr_centrality_measure_cache_.GetCachedResult(result, tdeps)) {
      SmartPtr<const Vector> compl_x_L = curr_compl_x_L();
      SmartPtr<const Vector> compl_x_U = curr_compl_x_U();
      SmartPtr<const Vector> compl_s_L = curr_compl_s_L();
      SmartPtr<const Vector> compl_s_U = curr_compl_s_U();

      result = CalcCentralityMeasure(*compl_x_L, *compl_x_U,
                                     *compl_s_L, *compl_s_U);

      curr_centrality_measure_cache_.AddCachedResult(result, tdeps);
    }
    return result;
  }

  Number
  IpoptCalculatedQuantities::curr_nlp_error()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_nlp_error()",
                   dbg_verbosity);
    DBG_ASSERT(initialize_called_);
    Number result;

    SmartPtr<const Vector> x = ip_data_->curr_x();
    SmartPtr<const Vector> s = ip_data_->curr_s();
    SmartPtr<const Vector> y_c = ip_data_->curr_y_c();
    SmartPtr<const Vector> y_d = ip_data_->curr_y_d();
    SmartPtr<const Vector> z_L = ip_data_->curr_z_L();
    SmartPtr<const Vector> z_U = ip_data_->curr_z_U();
    SmartPtr<const Vector> v_L = ip_data_->curr_v_L();
    SmartPtr<const Vector> v_U = ip_data_->curr_v_U();

    std::vector<const TaggedObject*> tdeps;
    tdeps.push_back(GetRawPtr(x));
    tdeps.push_back(GetRawPtr(s));
    tdeps.push_back(GetRawPtr(y_c));
    tdeps.push_back(GetRawPtr(y_c));
    tdeps.push_back(GetRawPtr(z_L));
    tdeps.push_back(GetRawPtr(z_U));
    tdeps.push_back(GetRawPtr(v_L));
    tdeps.push_back(GetRawPtr(v_U));

    if (!curr_nlp_error_cache_.GetCachedResult(result, tdeps)) {
      Number s_d = 0;
      Number s_c = 0;
      ComputeOptimalityErrorScaling(*ip_data_->curr_y_c(), *ip_data_->curr_y_d(),
                                    *ip_data_->curr_z_L(), *ip_data_->curr_z_U(),
                                    *ip_data_->curr_v_L(), *ip_data_->curr_v_U(),
                                    s_max_,
                                    s_d, s_c);
      DBG_PRINT((1, "s_d = %lf, s_c = %lf\n", s_d, s_c));

      // Primal infeasibility
      DBG_PRINT((1, "curr_dual_infeasibility(NORM_MAX) = %8.2e\n",
		 curr_dual_infeasibility(NORM_MAX)));
      result = curr_dual_infeasibility(NORM_MAX)/s_d;
      // Dual infeasibility
      DBG_PRINT((1, "curr_primal_infeasibility(NORM_MAX) = %8.2e\n",
		 curr_primal_infeasibility(NORM_MAX)));
      result = Max(result, curr_primal_infeasibility(NORM_MAX));
      // Complementarity
      DBG_PRINT((1, "curr_complementarity(0., NORM_MAX) = %8.2e\n",
		 curr_complementarity(0., NORM_MAX)));
      result = Max(result, curr_complementarity(0., NORM_MAX)/s_c);

      curr_nlp_error_cache_.AddCachedResult(result, tdeps);
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::curr_barrier_error()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_barrier_error()",
                   dbg_verbosity);
    DBG_ASSERT(initialize_called_);
    Number result;

    SmartPtr<const Vector> x = ip_data_->curr_x();
    SmartPtr<const Vector> s = ip_data_->curr_s();
    SmartPtr<const Vector> y_c = ip_data_->curr_y_c();
    SmartPtr<const Vector> y_d = ip_data_->curr_y_d();
    SmartPtr<const Vector> z_L = ip_data_->curr_z_L();
    SmartPtr<const Vector> z_U = ip_data_->curr_z_U();
    SmartPtr<const Vector> v_L = ip_data_->curr_v_L();
    SmartPtr<const Vector> v_U = ip_data_->curr_v_U();
    Number mu = ip_data_->curr_mu();

    std::vector<const TaggedObject*> tdeps;
    tdeps.push_back(GetRawPtr(x));
    tdeps.push_back(GetRawPtr(s));
    tdeps.push_back(GetRawPtr(y_c));
    tdeps.push_back(GetRawPtr(y_c));
    tdeps.push_back(GetRawPtr(z_L));
    tdeps.push_back(GetRawPtr(z_U));
    tdeps.push_back(GetRawPtr(v_L));
    tdeps.push_back(GetRawPtr(v_U));
    std::vector<Number> sdeps;
    sdeps.push_back(mu);

    if (!curr_barrier_error_cache_.GetCachedResult(result, tdeps, sdeps)) {
      Number s_d = 0;
      Number s_c = 0;
      ComputeOptimalityErrorScaling(*ip_data_->curr_y_c(), *ip_data_->curr_y_d(),
                                    *ip_data_->curr_z_L(), *ip_data_->curr_z_U(),
                                    *ip_data_->curr_v_L(), *ip_data_->curr_v_U(),
                                    s_max_,
                                    s_d, s_c);
      DBG_PRINT((1, "s_d = %lf, s_c = %lf\n", s_d, s_c));

      // Primal infeasibility
      result = curr_dual_infeasibility(NORM_MAX)/s_d;
      // Dual infeasibility
      result = Max(result, curr_primal_infeasibility(NORM_MAX));
      // Complementarity
      result = Max(result, curr_complementarity(mu, NORM_MAX)/s_c);

      curr_barrier_error_cache_.AddCachedResult(result, tdeps);
    }

    return result;
  }

  /*
  Number
  IpoptCalculatedQuantities::curr_primal_dual_error()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_primal_dual_error()",
                   dbg_verbosity);
    DBG_ASSERT(initialize_called_);
    Number result;

    SmartPtr<const Vector> x = ip_data_->curr_x();
    SmartPtr<const Vector> s = ip_data_->curr_s();
    SmartPtr<const Vector> y_c = ip_data_->curr_y_c();
    SmartPtr<const Vector> y_d = ip_data_->curr_y_d();
    SmartPtr<const Vector> z_L = ip_data_->curr_z_L();
    SmartPtr<const Vector> z_U = ip_data_->curr_z_U();
    SmartPtr<const Vector> v_L = ip_data_->curr_v_L();
    SmartPtr<const Vector> v_U = ip_data_->curr_v_U();

    std::vector<const TaggedObject*> tdeps;
    tdeps.push_back(GetRawPtr(x));
    tdeps.push_back(GetRawPtr(s));
    tdeps.push_back(GetRawPtr(y_c));
    tdeps.push_back(GetRawPtr(y_c));
    tdeps.push_back(GetRawPtr(z_L));
    tdeps.push_back(GetRawPtr(z_U));
    tdeps.push_back(GetRawPtr(v_L));
    tdeps.push_back(GetRawPtr(v_U));

    if (!curr_primal_dual_error_cache_.GetCachedResult(result, tdeps)) {
      Number s_d = 0;
      Number s_c = 0;
      ComputeOptimalityErrorScaling(*ip_data_->curr_y_c(), *ip_data_->curr_y_d(),
                                    *ip_data_->curr_z_L(), *ip_data_->curr_z_U(),
                                    *ip_data_->curr_v_L(), *ip_data_->curr_v_U(),
                                    s_max_,
                                    s_d, s_c);
      DBG_PRINT((1, "s_d = %lf, s_c = %lf\n", s_d, s_c));
      DBG_PRINT_VECTOR(2,"curr_grad_lag_x", *curr_grad_lag_x());
      result = curr_grad_lag_x()->Amax()/s_d;
      DBG_PRINT((1, "result = %lf\n", result));
      DBG_PRINT_VECTOR(2,"curr_grad_lag_s", *curr_grad_lag_s());
      result = Max(result, curr_grad_lag_s()->Amax()/s_d);
      DBG_PRINT((1, "result = %lf\n", result));

      DBG_PRINT_VECTOR(2,"curr_c", *curr_c());
      result = Max(result, curr_c()->Amax());
      DBG_PRINT((1, "result = %lf\n", result));
      DBG_PRINT_VECTOR(2,"curr_d_minus_s", *curr_d_minus_s());
      result = Max(result, curr_d_minus_s()->Amax());
      DBG_PRINT((1, "result = %lf\n", result));

      DBG_PRINT_VECTOR(2,"curr_compl_x_l", *curr_compl_x_L());
      result = Max(result, curr_compl_x_L()->Amax()/s_c);
      DBG_PRINT((1, "result = %lf\n", result));
      DBG_PRINT_VECTOR(2,"curr_compl_x_U", *curr_compl_x_U());
      result = Max(result, curr_compl_x_U()->Amax()/s_c);
      DBG_PRINT((1, "result = %lf\n", result));
      DBG_PRINT_VECTOR(2,"curr_compl_s_l", *curr_compl_s_L());
      result = Max(result, curr_compl_s_L()->Amax()/s_c);
      DBG_PRINT((1, "result = %lf\n", result));
      DBG_PRINT_VECTOR(2,"curr_compl_s_U", *curr_compl_s_U());
      result = Max(result, curr_compl_s_U()->Amax()/s_c);
      DBG_PRINT((1, "result = %lf\n", result));
      curr_primal_dual_error_cache_.AddCachedResult(result, tdeps);
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::curr_relaxed_primal_dual_error()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_relaxed_primal_dual_error()",
                   dbg_verbosity);
    DBG_ASSERT(initialize_called_);
    Number result;

    SmartPtr<const Vector> x = ip_data_->curr_x();
    SmartPtr<const Vector> s = ip_data_->curr_s();
    SmartPtr<const Vector> y_c = ip_data_->curr_y_c();
    SmartPtr<const Vector> y_d = ip_data_->curr_y_d();
    SmartPtr<const Vector> z_L = ip_data_->curr_z_L();
    SmartPtr<const Vector> z_U = ip_data_->curr_z_U();
    SmartPtr<const Vector> v_L = ip_data_->curr_v_L();
    SmartPtr<const Vector> v_U = ip_data_->curr_v_U();

    std::vector<const TaggedObject*> tdeps;
    tdeps.push_back(GetRawPtr(x));
    tdeps.push_back(GetRawPtr(s));
    tdeps.push_back(GetRawPtr(y_c));
    tdeps.push_back(GetRawPtr(y_c));
    tdeps.push_back(GetRawPtr(z_L));
    tdeps.push_back(GetRawPtr(z_U));
    tdeps.push_back(GetRawPtr(v_L));
    tdeps.push_back(GetRawPtr(v_U));

    Number mu = ip_data_->curr_mu();
    std::vector<Number> sdeps;
    sdeps.push_back(mu);

    if (!curr_relaxed_primal_dual_error_cache_.GetCachedResult(result, tdeps, sdeps)) {
      Number s_d = 0;
      Number s_c = 0;
      ComputeOptimalityErrorScaling(*ip_data_->curr_y_c(), *ip_data_->curr_y_d(),
                                    *ip_data_->curr_z_L(), *ip_data_->curr_z_U(),
                                    *ip_data_->curr_v_L(), *ip_data_->curr_v_U(),
                                    s_max_,
                                    s_d, s_c);

      result = curr_grad_lag_x()->Amax()/s_d;
      result = Max(result, curr_grad_lag_s()->Amax()/s_d);
      result = Max(result, curr_c()->Amax());
      result = Max(result, curr_d_minus_s()->Amax());
      result = Max(result, curr_relaxed_compl_x_L()->Amax()/s_c);
      result = Max(result, curr_relaxed_compl_x_U()->Amax()/s_c);
      result = Max(result, curr_relaxed_compl_s_L()->Amax()/s_c);
      result = Max(result, curr_relaxed_compl_s_U()->Amax()/s_c);
      curr_relaxed_primal_dual_error_cache_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;
  }
  */

  ///////////////////////////////////////////////////////////////////////////
  //                Fraction-to-the-boundary step sizes                    //
  ///////////////////////////////////////////////////////////////////////////

  Number
  IpoptCalculatedQuantities::CalcFracToZeroBound(const Vector& x,
      const Vector& delta,
      Number tau)
  {
    DBG_START_METH("IpoptCalculatedQuantities::CalcFracToZeroBound",
                   dbg_verbosity);
    DBG_ASSERT(x.Dim() == delta.Dim());
    if (x.Dim() == 0 && delta.Dim() == 0) {
      return 1.0;
    }

    SmartPtr<Vector> inv_alpha_bar = x.MakeNew();
    inv_alpha_bar->Copy(delta);
    inv_alpha_bar->Scal(-1.0/tau);
    inv_alpha_bar->ElementWiseDivide(x);
    DBG_PRINT_VECTOR(2, "x", x);
    DBG_PRINT_VECTOR(2, "inv_alpha_bar", *inv_alpha_bar);

    Number alpha = inv_alpha_bar->Max();
    if (alpha > 0) {
      alpha = Min(1.0/alpha, 1.0);
    }
    else {
      alpha = 1.0;
    }

    return alpha;
  }

  Number
  IpoptCalculatedQuantities::CalcFracToBound(const Vector& slack_L,
      const Matrix& P_L,
      const Vector& slack_U,
      const Matrix& P_U,
      const Vector& delta,
      Number tau)
  {
    DBG_START_METH("IpoptCalculatedQuantities::CalcFracToBound",
                   dbg_verbosity);

    Number alpha_L = 1.0;
    Number alpha_U = 1.0;
    if (slack_L.Dim() > 0) {
      SmartPtr<Vector> compressed_delta = slack_L.MakeNew();
      P_L.TransMultVector(1.0, delta, 0.0, *compressed_delta);

      SmartPtr<Vector> inv_alpha_bar = slack_L.MakeNew();
      inv_alpha_bar->Copy(*compressed_delta);
      inv_alpha_bar->Scal(-1.0/tau);
      inv_alpha_bar->ElementWiseDivide(slack_L);

      alpha_L = inv_alpha_bar->Max();
      if (alpha_L > 0) {
        alpha_L = Min((1.0/alpha_L), 1.0);
      }
      else {
        alpha_L = 1.0;
      }
    }

    if (slack_U.Dim() > 0) {
      SmartPtr<Vector> compressed_delta = slack_U.MakeNew();
      P_U.TransMultVector(1.0, delta, 0.0, *compressed_delta);

      SmartPtr<Vector> inv_alpha_bar = slack_U.MakeNew();
      inv_alpha_bar->Copy(*compressed_delta);
      inv_alpha_bar->Scal(1.0/tau);
      inv_alpha_bar->ElementWiseDivide(slack_U);

      alpha_U = inv_alpha_bar->Max();
      if (alpha_U > 0) {
        alpha_U = Min((1.0/alpha_U), 1.0);
      }
      else {
        alpha_U = 1.0;
      }
    }

    DBG_PRINT((1,"alpha_L = %lf, alpha_U = %lf\n", alpha_L, alpha_U));
    DBG_ASSERT(alpha_L >= 0.0 && alpha_L <= 1.0
               && alpha_U >=0.0 && alpha_U <= 1.0);

    return Min(alpha_L, alpha_U);
  }

  Number
  IpoptCalculatedQuantities::primal_frac_to_the_bound(Number tau,
      const Vector& delta_x,
      const Vector& delta_s)
  {
    DBG_START_METH("IpoptCalculatedQuantities::primal_frac_to_the_bound",
                   dbg_verbosity);
    Number result;
    SmartPtr<const Vector> x = ip_data_->curr_x();
    SmartPtr<const Vector> s = ip_data_->curr_s();
    std::vector<const TaggedObject*> tdeps;
    tdeps.push_back(GetRawPtr(x));
    tdeps.push_back(GetRawPtr(s));
    tdeps.push_back(&delta_x);
    tdeps.push_back(&delta_s);

    std::vector<Number> sdeps;
    sdeps.push_back(tau);

    if (!primal_frac_to_the_bound_cache_.GetCachedResult(result, tdeps, sdeps)) {
      result = Min(CalcFracToBound(*curr_slack_x_L(), *ip_nlp_->Px_L(),
                                   *curr_slack_x_U(), *ip_nlp_->Px_U(),
                                   delta_x, tau),
                   CalcFracToBound(*curr_slack_s_L(), *ip_nlp_->Pd_L(),
                                   *curr_slack_s_U(), *ip_nlp_->Pd_U(),
                                   delta_s, tau));

      primal_frac_to_the_bound_cache_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::curr_primal_frac_to_the_bound(Number tau)
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_primal_frac_to_the_bound()",
                   dbg_verbosity);
    return primal_frac_to_the_bound(tau, *ip_data_->delta_x(),
                                    *ip_data_->delta_s());
  }

  Number
  IpoptCalculatedQuantities::dual_frac_to_the_bound(Number tau,
      const Vector& delta_z_L,
      const Vector& delta_z_U,
      const Vector& delta_v_L,
      const Vector& delta_v_U)
  {
    DBG_START_METH("IpoptCalculatedQuantities::dual_frac_to_the_bound",
                   dbg_verbosity);
    Number result;

    SmartPtr<const Vector> z_L = ip_data_->curr_z_L();
    SmartPtr<const Vector> z_U = ip_data_->curr_z_U();
    SmartPtr<const Vector> v_L = ip_data_->curr_v_L();
    SmartPtr<const Vector> v_U = ip_data_->curr_v_U();
    std::vector<const TaggedObject*> tdeps;
    tdeps.push_back(GetRawPtr(z_L));
    tdeps.push_back(GetRawPtr(z_U));
    tdeps.push_back(GetRawPtr(v_L));
    tdeps.push_back(GetRawPtr(v_U));
    tdeps.push_back(&delta_z_L);
    tdeps.push_back(&delta_z_U);
    tdeps.push_back(&delta_v_L);
    tdeps.push_back(&delta_v_U);

    std::vector<Number> sdeps;
    sdeps.push_back(tau);

    if (!dual_frac_to_the_bound_cache_.GetCachedResult(result, tdeps, sdeps)) {
      result = CalcFracToZeroBound(*z_L, delta_z_L, tau);
      result = Min(result, CalcFracToZeroBound(*z_U, delta_z_U, tau));
      result = Min(result, CalcFracToZeroBound(*v_L, delta_v_L, tau));
      result = Min(result, CalcFracToZeroBound(*v_U, delta_v_U, tau));

      dual_frac_to_the_bound_cache_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::curr_dual_frac_to_the_bound(Number tau)
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_dual_frac_to_the_bound()",
                   dbg_verbosity);
    return dual_frac_to_the_bound(tau, *ip_data_->delta_z_L(),
                                  *ip_data_->delta_z_U(),
                                  *ip_data_->delta_v_L(),
                                  *ip_data_->delta_v_U());
  }

  Number
  IpoptCalculatedQuantities::slack_frac_to_the_bound(Number tau,
      const Vector& delta_x_L,
      const Vector& delta_x_U,
      const Vector& delta_s_L,
      const Vector& delta_s_U)
  {
    DBG_START_METH("IpoptCalculatedQuantities::slack_frac_to_the_bound",
                   dbg_verbosity);
    Number result;

    SmartPtr<const Vector> x_L = curr_slack_x_L();
    SmartPtr<const Vector> x_U = curr_slack_x_U();
    SmartPtr<const Vector> s_L = curr_slack_s_L();
    SmartPtr<const Vector> s_U = curr_slack_s_U();
    std::vector<const TaggedObject*> tdeps;
    tdeps.push_back(GetRawPtr(x_L));
    tdeps.push_back(GetRawPtr(x_U));
    tdeps.push_back(GetRawPtr(s_L));
    tdeps.push_back(GetRawPtr(s_U));
    tdeps.push_back(&delta_x_L);
    tdeps.push_back(&delta_x_U);
    tdeps.push_back(&delta_s_L);
    tdeps.push_back(&delta_s_U);

    std::vector<Number> sdeps;
    sdeps.push_back(tau);

    if (!slack_frac_to_the_bound_cache_.GetCachedResult(result, tdeps, sdeps)) {
      result = CalcFracToZeroBound(*x_L, delta_x_L, tau);
      result = Min(result, CalcFracToZeroBound(*x_U, delta_x_U, tau));
      result = Min(result, CalcFracToZeroBound(*s_L, delta_s_L, tau));
      result = Min(result, CalcFracToZeroBound(*s_U, delta_s_U, tau));

      slack_frac_to_the_bound_cache_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;
  }

  ///////////////////////////////////////////////////////////////////////////
  //                             Sigma Matrices                            //
  ///////////////////////////////////////////////////////////////////////////

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_sigma_x()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_sigma_x()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->curr_x();
    SmartPtr<const Vector> z_L = ip_data_->curr_z_L();
    SmartPtr<const Vector> z_U = ip_data_->curr_z_U();

    if (!curr_sigma_x_cache_.GetCachedResult3Dep(result, *x, *z_L, *z_U)) {
      SmartPtr<Vector> sigma = x->MakeNew();

      SmartPtr<const Vector> slack = curr_slack_x_L();
      SmartPtr<Vector> tmp = slack->MakeNew();
      tmp->Copy(*z_L);
      DBG_PRINT_VECTOR(2,"z_L", *z_L);
      DBG_PRINT_VECTOR(2,"slack_x_L", *slack);
      tmp->ElementWiseDivide(*slack);
      DBG_PRINT_VECTOR(2,"z_L/slack_x_L", *tmp);
      ip_nlp_->Px_L()->MultVector(1.0, *tmp, 0.0, *sigma);

      slack =  curr_slack_x_U();
      tmp = slack->MakeNew();
      tmp->Copy(*z_U);
      tmp->ElementWiseDivide(*slack);
      ip_nlp_->Px_U()->MultVector(1.0, *tmp, 1.0, *sigma);

      DBG_PRINT_VECTOR(2,"sigma_x", *sigma);

      result = ConstPtr(sigma);
      curr_sigma_x_cache_.AddCachedResult3Dep(result, *x, *z_L, *z_U);
    }

    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_sigma_s()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_sigma_s()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> s = ip_data_->curr_s();
    SmartPtr<const Vector> v_L = ip_data_->curr_v_L();
    SmartPtr<const Vector> v_U = ip_data_->curr_v_U();

    if (!curr_sigma_s_cache_.GetCachedResult3Dep(result, *s, *v_L, *v_U)) {
      SmartPtr<Vector> sigma = s->MakeNew();

      SmartPtr<const Vector> slack = curr_slack_s_L();
      SmartPtr<Vector> tmp = slack->MakeNew();
      tmp->Copy(*v_L);
      tmp->ElementWiseDivide(*slack);
      ip_nlp_->Pd_L()->MultVector(1.0, *tmp, 0.0, *sigma);

      slack =  curr_slack_s_U();
      tmp = slack->MakeNew();
      tmp->Copy(*v_U);
      tmp->ElementWiseDivide(*slack);
      ip_nlp_->Pd_U()->MultVector(1.0, *tmp, 1.0, *sigma);

      result = ConstPtr(sigma);
      curr_sigma_s_cache_.AddCachedResult3Dep(result, *s, *v_L, *v_U);
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::curr_avrg_compl()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_avrg_compl()",
                   dbg_verbosity);

    Number result;

    SmartPtr<const Vector> x = ip_data_->curr_x();
    SmartPtr<const Vector> s = ip_data_->curr_s();
    SmartPtr<const Vector> z_L = ip_data_->curr_z_L();
    SmartPtr<const Vector> z_U = ip_data_->curr_z_U();
    SmartPtr<const Vector> v_L = ip_data_->curr_v_L();
    SmartPtr<const Vector> v_U = ip_data_->curr_v_U();

    std::vector<const TaggedObject*> tdeps;
    tdeps.push_back(GetRawPtr(x));
    tdeps.push_back(GetRawPtr(s));
    tdeps.push_back(GetRawPtr(z_L));
    tdeps.push_back(GetRawPtr(z_U));
    tdeps.push_back(GetRawPtr(v_L));
    tdeps.push_back(GetRawPtr(z_U));

    if (!curr_avrg_compl_cache_.GetCachedResult(result, tdeps)) {
      if (!trial_avrg_compl_cache_.GetCachedResult(result, tdeps)) {

        SmartPtr<const Vector> slack_x_L = curr_slack_x_L();
        SmartPtr<const Vector> slack_x_U = curr_slack_x_U();
        SmartPtr<const Vector> slack_s_L = curr_slack_s_L();
        SmartPtr<const Vector> slack_s_U = curr_slack_s_U();

        Index ncomps = z_L->Dim() + z_U->Dim() + v_L->Dim() + v_U->Dim();

        if (ncomps>0) {
          result = z_L->Dot(*slack_x_L);
          result += z_U->Dot(*slack_x_U);
          result += v_L->Dot(*slack_s_L);
          result += v_U->Dot(*slack_s_U);

          result /= (Number)ncomps;
        }
        else {
          result = 0.;
        }
      }

      curr_avrg_compl_cache_.AddCachedResult(result, tdeps);
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::trial_avrg_compl()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_avrg_compl()",
                   dbg_verbosity);

    Number result;

    SmartPtr<const Vector> x = ip_data_->trial_x();
    SmartPtr<const Vector> s = ip_data_->trial_s();
    SmartPtr<const Vector> z_L = ip_data_->trial_z_L();
    SmartPtr<const Vector> z_U = ip_data_->trial_z_U();
    SmartPtr<const Vector> v_L = ip_data_->trial_v_L();
    SmartPtr<const Vector> v_U = ip_data_->trial_v_U();

    std::vector<const TaggedObject*> tdeps;
    tdeps.push_back(GetRawPtr(x));
    tdeps.push_back(GetRawPtr(s));
    tdeps.push_back(GetRawPtr(z_L));
    tdeps.push_back(GetRawPtr(z_U));
    tdeps.push_back(GetRawPtr(v_L));
    tdeps.push_back(GetRawPtr(z_U));

    if (!trial_avrg_compl_cache_.GetCachedResult(result, tdeps)) {
      if (!curr_avrg_compl_cache_.GetCachedResult(result, tdeps)) {

        SmartPtr<const Vector> slack_x_L = trial_slack_x_L();
        SmartPtr<const Vector> slack_x_U = trial_slack_x_U();
        SmartPtr<const Vector> slack_s_L = trial_slack_s_L();
        SmartPtr<const Vector> slack_s_U = trial_slack_s_U();

        Index ncomps = z_L->Dim() + z_U->Dim() + v_L->Dim() + v_U->Dim();

        if (ncomps>0) {
          result = z_L->Dot(*slack_x_L);
          result += z_U->Dot(*slack_x_U);
          result += v_L->Dot(*slack_s_L);
          result += v_U->Dot(*slack_s_U);

          result /= (Number)ncomps;
        }
        else {
          result = 0.;
        }
      }

      trial_avrg_compl_cache_.AddCachedResult(result, tdeps);
    }

    return result;
  }

  void IpoptCalculatedQuantities::ComputeOptimalityErrorScaling(const Vector& y_c, const Vector& y_d,
      const Vector& z_L, const Vector& z_U,
      const Vector& v_L, const Vector& v_U,
      Number s_max,
      Number& s_d, Number& s_c)
  {
    DBG_ASSERT(initialize_called_);

    s_c = z_L.Asum() + z_U.Asum() + v_L.Asum() + v_U.Asum();
    Number n = (z_L.Dim() + z_U.Dim() + v_L.Dim() + v_U.Dim());
    if (n == 0) {
      s_c = 1.0;
    }
    else {
      s_c = s_c / n;
      s_c = Max(s_max, s_c)/s_max;
    }

    s_d = y_c.Asum() + y_d.Asum() + z_L.Asum() + z_U.Asum() + v_L.Asum() + v_U.Asum();
    n = (y_c.Dim() + y_d.Dim() + z_L.Dim() + z_U.Dim() + v_L.Dim() + v_U.Dim());
    if ( n == 0 ) {
      s_d = 1.0;
    }
    else {
      s_d = s_d / n;
      s_d = Max(s_max, s_d)/s_max;
    }
  }

  Number IpoptCalculatedQuantities::curr_gradBarrTDelta()
  {
    Number result;

    SmartPtr<const Vector> x = ip_data_->curr_x();
    SmartPtr<const Vector> s = ip_data_->curr_s();
    SmartPtr<const Vector> delta_x = ip_data_->delta_x();
    SmartPtr<const Vector> delta_s = ip_data_->delta_s();
    std::vector<const TaggedObject*> tdeps;
    tdeps.push_back(GetRawPtr(x));
    tdeps.push_back(GetRawPtr(s));
    tdeps.push_back(GetRawPtr(delta_x));
    tdeps.push_back(GetRawPtr(delta_s));

    if (!curr_gradBarrTDelta_cache_.GetCachedResult(result, tdeps)) {
      result = curr_grad_barrier_obj_x()->Dot(*delta_x) +
               curr_grad_barrier_obj_s()->Dot(*delta_s);

      curr_gradBarrTDelta_cache_.AddCachedResult(result, tdeps);
    }
    return result;
  }

  // ToDo we probably want to get rid of this:
  bool IpoptCalculatedQuantities::in_restoration_phase()
  {
    const RestoIpoptNLP* resto_nlp = dynamic_cast<const RestoIpoptNLP*>(GetRawPtr(ip_nlp_));
    return (resto_nlp!=NULL);
  }

} // namespace Ipopt
