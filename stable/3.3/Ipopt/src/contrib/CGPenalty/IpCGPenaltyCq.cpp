// Copyright (C) 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpIpoptCalculatedQuantities.cpp 988 2007-06-01 21:57:27Z andreasw $
//
// Authors:  Andreas Waechter             IBM    2007-06-04
//               derived from IpIpoptCalculatedQuantities.cpp

#include "IpCGPenaltyCq.hpp"
#include "IpCGPenaltyData.hpp"

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

      curr_penalty_function_cache_(2),
      trial_penalty_function_cache_(5),
      curr_direct_deriv_penalty_function_cache_(1),
      curr_fast_direct_deriv_penalty_function_cache_(1),
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
  {
  }

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
    Number penalty = ip_data_->CGPenData().curr_penalty();
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
    Number penalty = ip_data_->CGPenData().curr_penalty();
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
    SmartPtr<const Vector> dy_c = ip_data_->CGPenData().delta_cgpen()->y_c();
    SmartPtr<const Vector> dy_d = ip_data_->CGPenData().delta_cgpen()->y_d();
    SmartPtr<const Vector> dx = ip_data_->CGPenData().delta_cgpen()->x();
    SmartPtr<const Vector> ds = ip_data_->CGPenData().delta_cgpen()->s();
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
    Number penalty = ip_data_->CGPenData().curr_penalty();
    std::vector<Number> sdeps(2);
    sdeps[0] = mu;
    sdeps[1] = penalty;

    if (!curr_direct_deriv_penalty_function_cache_.GetCachedResult(result, tdeps, sdeps)) {
      SmartPtr<const Vector> c = ip_cq_->curr_c();
      SmartPtr<const Vector> d_minus_s = ip_cq_->curr_d_minus_s();
      result = ip_cq_->curr_grad_barrier_obj_x()->Dot(*dx) +
               ip_cq_->curr_grad_barrier_obj_s()->Dot(*ds);
      result -= penalty*ip_cq_->curr_primal_infeasibility(NORM_2);
      result += c->Dot(*y_c);
      result += c->Dot(*dy_c);
      result += d_minus_s->Dot(*y_d);
      result += d_minus_s->Dot(*dy_d);
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
    DBG_ASSERT(ip_data_->CGPenData().HaveCgPenDeltas());
    SmartPtr<const Vector> dy_c = ip_data_->CGPenData().delta_cgfast()->y_c();
    SmartPtr<const Vector> dy_d = ip_data_->CGPenData().delta_cgfast()->y_d();
    SmartPtr<const Vector> dx = ip_data_->CGPenData().delta_cgfast()->x();
    SmartPtr<const Vector> ds = ip_data_->CGPenData().delta_cgfast()->s();
    std::vector<const TaggedObject*> tdeps(6);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(s);
    tdeps[2] = GetRawPtr(dy_c);
    tdeps[3] = GetRawPtr(dy_d);
    tdeps[4] = GetRawPtr(dx);
    tdeps[5] = GetRawPtr(ds);
    Number mu = ip_data_->curr_mu();
    Number penalty = ip_data_->CGPenData().curr_penalty();
    std::vector<Number> sdeps(2);
    sdeps[0] = mu;
    sdeps[1] = penalty;

    if (!curr_fast_direct_deriv_penalty_function_cache_.GetCachedResult(result, tdeps, sdeps)) {
      SmartPtr<const Vector> c = ip_cq_->curr_c();
      SmartPtr<const Vector> d_minus_s = ip_cq_->curr_d_minus_s();
      result = ip_cq_->curr_grad_barrier_obj_x()->Dot(*dx) +
               ip_cq_->curr_grad_barrier_obj_s()->Dot(*ds);
      result -= penalty*ip_cq_->curr_primal_infeasibility(NORM_2);
      result += c->Dot(*dy_c);
      result += d_minus_s->Dot(*dy_d);
      curr_fast_direct_deriv_penalty_function_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    return result;
  }

  Number CGPenaltyCq::curr_cg_pert_fact()
  {
    DBG_START_METH("CGPenaltyCq::curr_cg_pert_fact()",
                   dbg_verbosity);

    Number result;

    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> s = ip_data_->curr()->s();
    std::vector<const TaggedObject*> tdeps(2);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(s);
    Number penalty = ip_data_->CGPenData().curr_penalty();
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

} // namespace Ipopt
