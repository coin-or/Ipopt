// Copyright (C) 2007, 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter             IBM    2007-06-04
//               derived from IpIpoptCalculatedQuantities.cpp

#include "IpCGPenaltyCq.hpp"
#include "IpCGPenaltyData.hpp"
#include "IpTripletHelper.hpp"

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
# endif
#endif

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  CGPenaltyCq::CGPenaltyCq(IpoptNLP* ip_nlp,
                           IpoptData* ip_data,
                           IpoptCalculatedQuantities* ip_cq)
      :
      ip_nlp_(ip_nlp),
      ip_data_(ip_data),
      ip_cq_(ip_cq),

      curr_fast_direct_deriv_penalty_function_cache_(1),
      curr_jac_cd_norm_cache_(1),
      curr_scaled_y_Amax_cache_(1),
      curr_added_y_nrm2_cache_(1),
      curr_penalty_function_cache_(2),
      trial_penalty_function_cache_(5),
      curr_direct_deriv_penalty_function_cache_(1),
      curr_cg_pert_fact_cache_(1),

      initialize_called_(false)
  {
    DBG_START_METH("CGPenaltyCq::CGPenaltyCq", dbg_verbosity);
    DBG_ASSERT(ip_nlp_);
    DBG_ASSERT(ip_data_);
    DBG_ASSERT(ip_cq_);
  }

  CGPenaltyCq::~CGPenaltyCq()
  {}

  void CGPenaltyCq::RegisterOptions(const SmartPtr<RegisteredOptions>& roptions)
  {}

  bool CGPenaltyCq::Initialize(const Journalist& jnlst,
                               const OptionsList& options,
                               const std::string& prefix)
  {
    initialize_called_ = true;
    return true;
  }

  //////////////////////////////////////////////////////
  //   Methods for the Chen-Goldfarb penalty function //
  //////////////////////////////////////////////////////

  Number
  CGPenaltyCq::curr_jac_cd_norm(Index nrm_type)
  {
    DBG_START_METH("CGPenaltyCq::curr_jac_cd_norm()",
                   dbg_verbosity);

    Number result;
    SmartPtr<const Matrix> jac_c = ip_cq_->curr_jac_c();
    Index nnz = TripletHelper::GetNumberEntries(*jac_c);
    Number* values = new Number[nnz];
    TripletHelper::FillValues(nnz, *jac_c, values);
    Index count = 1;
    result = 0.;
    for (Index i=1; i<nnz; i++) {
      if (nrm_type == 3) {
        result = Max(result, fabs(values[i]));
      }
      if (nrm_type == 1) {
        result += fabs(values[i]);
        count ++;
      }
    }
    delete [] values;
    SmartPtr<const Matrix> jac_d = ip_cq_->curr_jac_d();
    nnz = TripletHelper::GetNumberEntries(*jac_d);
    values = new Number[nnz];
    TripletHelper::FillValues(nnz, *jac_d, values);
    for (Index i=1; i<nnz; i++) {
      if (nrm_type == 3) {
        result = Max(result, fabs(values[i]));
      }
      if (nrm_type == 1) {
        result += fabs(values[i]);
        count ++;
      }
    }
    delete [] values;
    if (nrm_type == 1) {
      result = result/count;
    }
    return result;
  }

  Number
  CGPenaltyCq::curr_penalty_function()
  {
    DBG_START_METH("CGPenaltyCq::curr_penalty_function()",
                   dbg_verbosity);
    Number result;
    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> s = ip_data_->curr()->s();
    std::vector<const TaggedObject*> tdeps(2);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(s);
    Number mu = ip_data_->curr_mu();
    Number penalty = CGPenData().curr_penalty();
    std::vector<Number> sdeps(2);
    sdeps[0] = mu;
    sdeps[1] = penalty;
    if (!curr_penalty_function_cache_.GetCachedResult(result, tdeps, sdeps)) {
      if (!trial_penalty_function_cache_.GetCachedResult(result, tdeps, sdeps)) {
        result = ip_cq_->curr_barrier_obj();
        result += penalty*ip_cq_->curr_primal_infeasibility(NORM_2);
      }
      curr_penalty_function_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    DBG_ASSERT(IsFiniteNumber(result));
    return result;
  }

  Number
  CGPenaltyCq::trial_penalty_function()
  {
    DBG_START_METH("CGPenaltyCq::trial_penalty_function()",
                   dbg_verbosity);
    Number result;
    SmartPtr<const Vector> x = ip_data_->trial()->x();
    SmartPtr<const Vector> s = ip_data_->trial()->s();
    std::vector<const TaggedObject*> tdeps(2);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(s);
    Number mu = ip_data_->curr_mu();
    Number penalty = CGPenData().curr_penalty();
    std::vector<Number> sdeps(2);
    sdeps[0] = mu;
    sdeps[1] = penalty;
    if (!trial_penalty_function_cache_.GetCachedResult(result, tdeps, sdeps)) {
      if (!curr_penalty_function_cache_.GetCachedResult(result, tdeps, sdeps)) {
        result = ip_cq_->trial_barrier_obj();
        result += penalty*ip_cq_->trial_primal_infeasibility(NORM_2);
      }
      trial_penalty_function_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    DBG_ASSERT(IsFiniteNumber(result));
    return result;
  }



  Number CGPenaltyCq::curr_direct_deriv_penalty_function()
  {
    DBG_START_METH("CGPenaltyCq::curr_direct_deriv_penalty_function()",
                   dbg_verbosity);

    Number result;
    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> s = ip_data_->curr()->s();
    SmartPtr<const Vector> y_c = ip_data_->curr()->y_c();
    SmartPtr<const Vector> y_d = ip_data_->curr()->y_d();
    SmartPtr<const Vector> dy_c = CGPenData().delta_cgpen()->y_c();
    SmartPtr<const Vector> dy_d = CGPenData().delta_cgpen()->y_d();
    SmartPtr<const Vector> dx = CGPenData().delta_cgpen()->x();
    SmartPtr<const Vector> ds = CGPenData().delta_cgpen()->s();
    std::vector<const TaggedObject*> tdeps(8);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(s);
    tdeps[2] = GetRawPtr(y_c);
    tdeps[3] = GetRawPtr(y_d);
    tdeps[4] = GetRawPtr(dy_c);
    tdeps[5] = GetRawPtr(dy_d);
    tdeps[6] = GetRawPtr(dx);
    tdeps[7] = GetRawPtr(ds);
    Number mu = ip_data_->curr_mu();
    Number penalty = CGPenData().curr_penalty();
    std::vector<Number> sdeps(2);
    sdeps[0] = mu;
    sdeps[1] = penalty;
    if (!curr_direct_deriv_penalty_function_cache_.GetCachedResult(result, tdeps, sdeps)) {
      result = ip_cq_->curr_grad_barrier_obj_x()->Dot(*dx) +
               ip_cq_->curr_grad_barrier_obj_s()->Dot(*ds);
      Number curr_inf = ip_cq_->curr_primal_infeasibility(NORM_2);
      result -= penalty*curr_inf;
      if (curr_inf != 0.) {
        Number fac = penalty*CGPenData().CurrPenaltyPert()/curr_inf;
        SmartPtr<const Vector> c = ip_cq_->curr_c();
        SmartPtr<const Vector> d_minus_s = ip_cq_->curr_d_minus_s();
        Number result1 = c->Dot(*y_c);
        result1 += c->Dot(*dy_c);
        result1 += d_minus_s->Dot(*y_d);
        result1 += d_minus_s->Dot(*dy_d);
        result1 *= fac;
        result += result1;
      }
      curr_direct_deriv_penalty_function_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    return result;
  }

  Number CGPenaltyCq::curr_fast_direct_deriv_penalty_function()
  {
    DBG_START_METH("CGPenaltyCq::curr_fast_direct_deriv_penalty_function()",
                   dbg_verbosity);

    Number result;
    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> s = ip_data_->curr()->s();
    DBG_ASSERT(CGPenData().HaveCgPenDeltas());
    SmartPtr<const Vector> dy_c = CGPenData().delta_cgfast()->y_c();
    SmartPtr<const Vector> dy_d = CGPenData().delta_cgfast()->y_d();
    SmartPtr<const Vector> dx = CGPenData().delta_cgfast()->x();
    SmartPtr<const Vector> ds = CGPenData().delta_cgfast()->s();
    std::vector<const TaggedObject*> tdeps(6);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(s);
    tdeps[2] = GetRawPtr(dy_c);
    tdeps[3] = GetRawPtr(dy_d);
    tdeps[4] = GetRawPtr(dx);
    tdeps[5] = GetRawPtr(ds);
    Number mu = ip_data_->curr_mu();
    Number penalty = CGPenData().curr_penalty();
    std::vector<Number> sdeps(2);
    sdeps[0] = mu;
    sdeps[1] = penalty;
    if (!curr_fast_direct_deriv_penalty_function_cache_.GetCachedResult(result, tdeps, sdeps)) {
      result = ip_cq_->curr_grad_barrier_obj_x()->Dot(*dx) +
               ip_cq_->curr_grad_barrier_obj_s()->Dot(*ds);
      Number curr_inf = ip_cq_->curr_primal_infeasibility(NORM_2);
      result -= penalty*curr_inf;
      if (curr_inf != 0.) {
        Number fac = penalty*CGPenData().CurrPenaltyPert()/curr_inf;
        SmartPtr<const Vector> c = ip_cq_->curr_c();
        SmartPtr<const Vector> d_minus_s = ip_cq_->curr_d_minus_s();
        Number result1 = c->Dot(*dy_c);
        result1 += d_minus_s->Dot(*dy_d);
        result1 *= fac;
        result += result1;
      }
      curr_fast_direct_deriv_penalty_function_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    return result;
  }

  Number CGPenaltyCq::curr_cg_pert_fact()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_cg_pert_fact()",
                   dbg_verbosity);

    Number result;
    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> s = ip_data_->curr()->s();
    std::vector<const TaggedObject*> tdeps(2);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(s);
    Number penalty = CGPenData().curr_kkt_penalty();
    std::vector<Number> sdeps(1);
    sdeps[0] = penalty;
    DBG_ASSERT(penalty>0.);
    if (!curr_cg_pert_fact_cache_.GetCachedResult(result, tdeps, sdeps)) {
      Number eq_2nrm = ip_cq_->curr_primal_infeasibility(NORM_2);
      result = eq_2nrm/penalty;
      curr_cg_pert_fact_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    return result;
  }

  Number CGPenaltyCq::dT_times_barH_times_d()
  {
    DBG_START_METH("IpoptCalculatedQuantities::dT_times_barH_times_d()",
                   dbg_verbosity);
    Number result;
    SmartPtr<const Vector> d_x = CGPenData().delta_cgfast()->x();
    SmartPtr<const Vector> d_s = CGPenData().delta_cgfast()->s();
    SmartPtr<const Vector> y_c = ip_data_->curr()->y_c();
    SmartPtr<const Vector> y_d = ip_data_->curr()->y_d();
    SmartPtr<const Vector> dy_c = CGPenData().delta_cgfast()->y_c();
    SmartPtr<const Vector> dy_d = CGPenData().delta_cgfast()->y_d();
    SmartPtr<const Vector> c = ip_cq_->curr_c();
    SmartPtr<const Vector> d_minus_s = ip_cq_->curr_d_minus_s();
    Number deriv_barrier_dx = ip_cq_->curr_grad_barrier_obj_x()->Dot(*d_x);
    Number deriv_barrier_dx_ds = deriv_barrier_dx + ip_cq_->curr_grad_barrier_obj_s()->Dot(*d_s);
    Number penalty = CGPenData().curr_penalty();
    result = -y_c->Dot(*dy_c);
    result -= y_d->Dot(*dy_d);
    result *= curr_cg_pert_fact();
    result -= deriv_barrier_dx_ds;
    result += c->Dot(*y_c);
    result += d_minus_s->Dot(*y_d);
    result -= c->Dot(*dy_c);
    result -= d_minus_s->Dot(*dy_d);
    result += penalty*ip_cq_->curr_primal_infeasibility(NORM_2);

    return result;
  }


  Number CGPenaltyCq::compute_curr_cg_penalty(const Number pen_des_fact )
  {
    DBG_START_METH("CGPenaltyCq::compute_curr_cg_penalty()",
                   dbg_verbosity);

    SmartPtr<const Vector> d_x = ip_data_->delta()->x();
    SmartPtr<const Vector> d_s = ip_data_->delta()->s();
    SmartPtr<const Vector> y_c = ip_data_->curr()->y_c();
    SmartPtr<const Vector> y_d = ip_data_->curr()->y_d();
    SmartPtr<const Vector> dy_c = ip_data_->delta()->y_c();
    SmartPtr<const Vector> dy_d = ip_data_->delta()->y_d();

    // Compute delta barrier times (delta x, delta s)
    Number deriv_barrier_dx = ip_cq_->curr_grad_barrier_obj_x()->Dot(*d_x);
    Number deriv_barrier_dx_ds = deriv_barrier_dx + ip_cq_->curr_grad_barrier_obj_s()->Dot(*d_s);
    // Compute delta x times the damped Hessian times delta x
    SmartPtr<const Vector> tem_jac_cT_times_y_c =
      ip_cq_->curr_jac_cT_times_vec(*y_c);
    SmartPtr<const Vector> tem_jac_cT_times_dy_c =
      ip_cq_->curr_jac_cT_times_vec(*dy_c);
    SmartPtr<Vector> tem_jac_cT_times_y_c_plus_dy_c =
      tem_jac_cT_times_y_c->MakeNew();
    tem_jac_cT_times_y_c_plus_dy_c->AddTwoVectors(1.,*tem_jac_cT_times_y_c, 1.,
        *tem_jac_cT_times_dy_c, 0.);
    SmartPtr<const Vector> tem_jac_dT_times_y_d =
      ip_cq_->curr_jac_dT_times_vec(*y_d);
    SmartPtr<const Vector> tem_jac_dT_times_dy_d =
      ip_cq_->curr_jac_cT_times_vec(*dy_c);
    SmartPtr<Vector> tem_jac_dT_times_y_d_plus_dy_d =
      tem_jac_cT_times_y_c->MakeNew();
    tem_jac_dT_times_y_d_plus_dy_d->AddTwoVectors(1.,*tem_jac_dT_times_y_d, 1.,
        *tem_jac_dT_times_dy_d, 0.);
    Number d_xs_times_damped_Hessian_times_d_xs = -deriv_barrier_dx_ds;
    d_xs_times_damped_Hessian_times_d_xs +=
      -(tem_jac_cT_times_y_c_plus_dy_c->Dot(*d_x)
        +tem_jac_dT_times_y_d_plus_dy_d->Dot(*d_x)
        -y_d->Dot(*d_s)
        -dy_d->Dot(*d_s));
    Number dxs_nrm = pow(d_x->Nrm2(), 2.) + pow(d_s->Nrm2(), 2.);
    d_xs_times_damped_Hessian_times_d_xs = Max(1e-8*dxs_nrm,
                                           d_xs_times_damped_Hessian_times_d_xs);
    Number infeasibility = ip_cq_->curr_primal_infeasibility(NORM_2);
    Number penalty = 0.;
    if (infeasibility > 0.) {
      Number deriv_inf = 0.;
      Number fac = CGPenData().CurrPenaltyPert()/infeasibility;
      SmartPtr<const Vector> c = ip_cq_->curr_c();
      SmartPtr<const Vector> d_minus_s = ip_cq_->curr_d_minus_s();
      if (CGPenData().HaveCgFastDeltas()) {
        SmartPtr<const Vector> fast_dy_c = CGPenData().delta_cgfast()->y_c();
        SmartPtr<const Vector> fast_dy_d = CGPenData().delta_cgfast()->y_d();
        deriv_inf += c->Dot(*fast_dy_c);
        deriv_inf += d_minus_s->Dot(*fast_dy_d);
        deriv_inf *= fac;
        deriv_inf -= infeasibility;
      }
      else {
        SmartPtr<const Vector> cgpen_dy_c = CGPenData().delta_cgpen()->y_c();
        SmartPtr<const Vector> cgpen_dy_d = CGPenData().delta_cgpen()->y_d();
        deriv_inf += c->Dot(*cgpen_dy_c);
        deriv_inf += c->Dot(*y_c);
        deriv_inf += d_minus_s->Dot(*cgpen_dy_d);
        deriv_inf += d_minus_s->Dot(*y_d);
        deriv_inf *= fac;
        deriv_inf -= infeasibility;
      }
      penalty = -(deriv_barrier_dx_ds + pen_des_fact*
                  d_xs_times_damped_Hessian_times_d_xs)/
                (deriv_inf + pen_des_fact*infeasibility);
    }

    return penalty;
  }


  Number CGPenaltyCq::compute_curr_cg_penalty_scale()
  {
    DBG_START_METH("CGPenaltyCq::compute_curr_cg_penalty_scale()",
                   dbg_verbosity);
    Number penalty;
    Number infeasibility = ip_cq_->curr_primal_infeasibility(NORM_2);
    if (!CGPenData().NeverTryPureNewton()) {
      penalty = Min(1e13, infeasibility*1e9);
    }
    else {
      Number reference = (curr_jac_cd_norm(1) +
                          ip_cq_->curr_primal_infeasibility(NORM_1)/
                          (ip_data_->curr()->y_c()->Dim()+
                           ip_data_->curr()->y_d()->Dim()))/2.;
      if (CGPenData().restor_iter() == ip_data_->iter_count() ||
          ip_data_->iter_count() == 0) {
        reference_infeasibility_ = Min(1., infeasibility);
      }
      Number i = CGPenData().restor_counter();
      Number fac = 4*1e-2*pow(1e1, i);
      //Number fac = 1e-2;
      penalty = Min(1e4,infeasibility) / (reference*fac*
                                          pow(reference_infeasibility_, 1));
    }

    return penalty;
  }

  Number CGPenaltyCq::curr_scaled_y_Amax()
  {
    DBG_START_METH("CGPenaltyCq::curr_scaled_y_Amax()",
                   dbg_verbosity);

    Number result;
    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> y_c = ip_data_->curr()->y_c();
    SmartPtr<const Vector> y_d = ip_data_->curr()->y_d();
    std::vector<const TaggedObject*> deps(3);
    deps[0] = GetRawPtr(x);
    deps[1] = GetRawPtr(y_c);
    deps[2] = GetRawPtr(y_d);
    if (!curr_scaled_y_Amax_cache_.GetCachedResult(result, deps)) {
      result = Max(y_c->Amax(), y_d->Amax());
      result /= Max(1., ip_cq_->curr_grad_f()->Amax());
      curr_scaled_y_Amax_cache_.AddCachedResult(result, deps);
    }
    return result;
  }

  Number CGPenaltyCq::curr_added_y_nrm2()
  {
    DBG_START_METH("CGPenaltyCq::curr_added_y_nrm2()",
                   dbg_verbosity);

    Number result;
    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> y_c = ip_data_->curr()->y_c();
    SmartPtr<const Vector> y_d = ip_data_->curr()->y_d();
    std::vector<const TaggedObject*> deps(3);
    deps[0] = GetRawPtr(x);
    deps[1] = GetRawPtr(y_c);
    deps[2] = GetRawPtr(y_d);
    if (!curr_added_y_nrm2_cache_.GetCachedResult(result, deps)) {
      SmartPtr<Vector> y_c_plus_dy_c = ip_data_->delta()->y_c()->MakeNew();
      SmartPtr<Vector> y_d_plus_dy_d = ip_data_->delta()->y_d()->MakeNew();
      y_c_plus_dy_c->AddTwoVectors(1.,*ip_data_->delta()->y_c(),
                                   1.,*ip_data_->curr()->y_c(), 0.);
      y_d_plus_dy_d->AddTwoVectors(1.,*ip_data_->delta()->y_d(),
                                   1.,*ip_data_->curr()->y_d(), 0.);
      result = sqrt(pow(y_c_plus_dy_c->Nrm2(),2)
                    + pow(y_d_plus_dy_d->Nrm2(),2) );
      curr_added_y_nrm2_cache_.AddCachedResult(result, deps);
    }
    return result;
  }

} // namespace Ipopt
