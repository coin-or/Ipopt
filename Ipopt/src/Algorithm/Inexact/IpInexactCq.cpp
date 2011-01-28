// Copyright (C) 2008, 2011 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter             IBM    2008-08-31

#include "IpInexactCq.hpp"
#include "IpIpoptNLP.hpp"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif


  InexactCq::InexactCq(IpoptNLP* ip_nlp,
                       IpoptData* ip_data,
                       IpoptCalculatedQuantities* ip_cq)
      :
      ip_nlp_(ip_nlp),
      ip_data_(ip_data),
      ip_cq_(ip_cq),

      curr_jac_cdT_times_curr_cdminuss_cache_(1),
      curr_scaling_slacks_cache_(1),
      curr_slack_scaled_d_minus_s_cache_(1),
      curr_scaled_Ac_norm_cache_(1),
      slack_scaled_norm_cache_(6),
      curr_W_times_vec_x_cache_(0), // ToDo: decide if we want this cached
      curr_W_times_vec_s_cache_(0),
      curr_Wu_x_cache_(1),
      curr_Wu_s_cache_(1),
      curr_uWu_cache_(1),
      curr_jac_times_normal_c_cache_(1),
      curr_jac_times_normal_d_cache_(1)
  {
    DBG_START_METH("InexactCq::InexactCq", dbg_verbosity);
    DBG_ASSERT(ip_nlp_);
    DBG_ASSERT(ip_data_);
    DBG_ASSERT(ip_cq_);
  }

  InexactCq::~InexactCq()
  {}

  void
  InexactCq::RegisterOptions(const SmartPtr<RegisteredOptions>& roptions)
  {
    roptions->AddLowerBoundedNumberOption(
      "slack_scale_max",
      "Upper bound on slack-based scaling parameters.",
      0.0, true, 1.,
      "");
  }

  bool
  InexactCq::Initialize(const Journalist& jnlst,
                        const OptionsList& options,
                        const std::string& prefix)
  {
    options.GetNumericValue("slack_scale_max", slack_scale_max_, prefix);
    return true;
  }

  SmartPtr<const Vector>
  InexactCq::curr_jac_cdT_times_curr_cdminuss()
  {
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> s = ip_data_->curr()->s();

    if (!curr_jac_cdT_times_curr_cdminuss_cache_.GetCachedResult2Dep(result, *x, *s)) {
      // c part
      SmartPtr<const Vector> curr_c = ip_cq_->curr_c();
      SmartPtr<const Vector> curr_jac_cT_times_curr_c =
        ip_cq_->curr_jac_cT_times_vec(*curr_c);

      // d minus s part
      SmartPtr<const Vector> curr_d_minus_s = ip_cq_->curr_d_minus_s();
      SmartPtr<const Vector> curr_jac_dT_times_curr_d_minus_s =
        ip_cq_->curr_jac_dT_times_vec(*curr_d_minus_s);

      // add them
      SmartPtr<Vector> tmp = curr_jac_cT_times_curr_c->MakeNew();
      tmp->AddTwoVectors(1., *curr_jac_cT_times_curr_c,
                         1., *curr_jac_dT_times_curr_d_minus_s, 0.);
      result = ConstPtr(tmp);
      curr_jac_cdT_times_curr_cdminuss_cache_.AddCachedResult2Dep(result, *x, *s);
    }

    return result;
  }

  SmartPtr<const Vector>
  InexactCq::curr_scaling_slacks()
  {
    DBG_START_METH("InexactCq::curr_scaling_slacks()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> s = ip_data_->curr()->s();

    if (!curr_scaling_slacks_cache_.GetCachedResult1Dep(result, *s)) {
      SmartPtr<Vector> tmp = s->MakeNew();

      // Lower bounds
      SmartPtr<const Matrix> Pd_L = ip_nlp_->Pd_L();
      SmartPtr<const Vector> curr_slack_s_L = ip_cq_->curr_slack_s_L();
      DBG_PRINT_MATRIX(1, "Pd_L", *Pd_L);
      DBG_PRINT_VECTOR(1, "curr_slack_s_L", *curr_slack_s_L);
      Pd_L->MultVector(1., *curr_slack_s_L, 0., *tmp);

      // Upper bounds
      SmartPtr<const Matrix> Pd_U = ip_nlp_->Pd_U();
      SmartPtr<const Vector> curr_slack_s_U = ip_cq_->curr_slack_s_U();
      DBG_PRINT_MATRIX(1, "Pd_U", *Pd_U);
      DBG_PRINT_VECTOR(1, "curr_slack_s_U", *curr_slack_s_U);
      Pd_U->MultVector(1., *curr_slack_s_U, 1., *tmp);

      SmartPtr<Vector> tmp2 = tmp->MakeNew();
      tmp2->Set(slack_scale_max_);
      tmp->ElementWiseMin(*tmp2);

      result = ConstPtr(tmp);
      curr_scaling_slacks_cache_.AddCachedResult1Dep(result, *s);
    }

    return result;
  }

  SmartPtr<const Vector>
  InexactCq::curr_slack_scaled_d_minus_s()
  {
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> s = ip_data_->curr()->s();

    if (!curr_slack_scaled_d_minus_s_cache_.GetCachedResult2Dep(result, *x, *s)) {
      SmartPtr<const Vector> d_minus_s = ip_cq_->curr_d_minus_s();
      SmartPtr<const Vector> scaling_slacks = curr_scaling_slacks();

      SmartPtr<Vector> tmp = d_minus_s->MakeNewCopy();
      tmp->ElementWiseMultiply(*scaling_slacks);
      result = ConstPtr(tmp);
      curr_slack_scaled_d_minus_s_cache_.AddCachedResult2Dep(result, *x, *s);
    }

    return result;
  }

  Number
  InexactCq::curr_scaled_Ac_norm()
  {
    DBG_START_METH("InexactCq::curr_scaled_Ac_norm",
                   dbg_verbosity);
    Number result;

    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> s = ip_data_->curr()->s();

    if (!curr_scaled_Ac_norm_cache_.GetCachedResult2Dep(result, *x, *s)) {
      SmartPtr<const Vector> jac_cdT_times_curr_cdminuss =
        curr_jac_cdT_times_curr_cdminuss();
      DBG_PRINT_VECTOR(2, "jac_cdT_times_curr_cdminuss",
                       *jac_cdT_times_curr_cdminuss);
      SmartPtr<const Vector> slack_scaled_d_minus_s =
        curr_slack_scaled_d_minus_s();
      DBG_PRINT_VECTOR(2, "slack_scaled_d_minus_s",
                       *slack_scaled_d_minus_s);
      result = ip_cq_->CalcNormOfType(NORM_2, *jac_cdT_times_curr_cdminuss,
                                      *slack_scaled_d_minus_s);

      curr_scaled_Ac_norm_cache_.AddCachedResult2Dep(result, *x, *s);
    }

    return result;
  }

  /** ||A||^2 */
  Number
  InexactCq::curr_scaled_A_norm2()
  {
    DBG_START_METH("InexactCq::curr_scaled_A_norm",
                   dbg_verbosity);
    Number result;

    //if (!curr_scaled_A_norm_cache_.GetCachedResult(...)) {
    result = 2;
    //curr_scaled_A_norm_cache_.AddCachedResult(...);
    //}

    return result;
  }

  Number
  InexactCq::slack_scaled_norm(const Vector& x, const Vector &s)
  {
    SmartPtr<const Vector> scaling_slacks = curr_scaling_slacks();

    Number result;
    if (!slack_scaled_norm_cache_.GetCachedResult3Dep(result, *scaling_slacks, x, s)) {
      SmartPtr<Vector> tmp = s.MakeNewCopy();
      tmp->ElementWiseDivide(*scaling_slacks);
      result = ip_cq_->CalcNormOfType(NORM_2, x, *tmp);
      slack_scaled_norm_cache_.AddCachedResult3Dep(result, *scaling_slacks, x, s);
    }

    return result;
  }

  SmartPtr<const Vector>
  InexactCq::curr_W_times_vec_x(const Vector& vec_x)
  {
    DBG_START_METH("InexactCq::curr_W_times_vec_x",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const SymMatrix> W = ip_data_->W();
    Number pd_pert_x;
    Number pd_pert_s;
    Number pd_pert_c;
    Number pd_pert_d;
    ip_data_->getPDPert(pd_pert_x, pd_pert_s, pd_pert_c, pd_pert_d);

    std::vector<const TaggedObject*> tdeps(2);
    tdeps[0] = &vec_x;
    tdeps[1] = GetRawPtr(W);
    std::vector<Number> sdeps(1);
    sdeps[0] = pd_pert_x;

    if (!curr_W_times_vec_x_cache_.GetCachedResult(result, tdeps, sdeps)) {
      DBG_PRINT_VECTOR(2, "vec_x", vec_x);
      DBG_PRINT_MATRIX(2, "W", *W);
      DBG_PRINT((2, "pd_pert_x = %e\n", pd_pert_x));
      SmartPtr<Vector> tmp = vec_x.MakeNewCopy();
      W->MultVector(1., vec_x, pd_pert_x, *tmp);
      result = ConstPtr(tmp);
      curr_W_times_vec_x_cache_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;
  }

  SmartPtr<const Vector>
  InexactCq::curr_W_times_vec_s(const Vector& vec_s)
  {
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> sigma_s = ip_cq_->curr_sigma_s();
    Number pd_pert_x;
    Number pd_pert_s;
    Number pd_pert_c;
    Number pd_pert_d;
    ip_data_->getPDPert(pd_pert_x, pd_pert_s, pd_pert_c, pd_pert_d);

    std::vector<const TaggedObject*> tdeps(2);
    tdeps[0] = &vec_s;
    tdeps[1] = GetRawPtr(sigma_s);
    std::vector<Number> sdeps(1);
    sdeps[0] = pd_pert_s;

    if (!curr_W_times_vec_s_cache_.GetCachedResult(result, tdeps, sdeps)) {
      SmartPtr<Vector> tmp = vec_s.MakeNewCopy();
      tmp->ElementWiseMultiply(*sigma_s);
      if (pd_pert_s>0.) {
        tmp->AddOneVector(pd_pert_s, vec_s, 1.);
      }
      result = ConstPtr(tmp);
      curr_W_times_vec_s_cache_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;
  }

  SmartPtr<const Vector>
  InexactCq::curr_Wu_x()
  {
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> tangential_x = InexData().tangential_x();
    SmartPtr<const SymMatrix> W = ip_data_->W();
    Number pd_pert_x;
    Number pd_pert_s;
    Number pd_pert_c;
    Number pd_pert_d;
    ip_data_->getPDPert(pd_pert_x, pd_pert_s, pd_pert_c, pd_pert_d);

    std::vector<const TaggedObject*> tdeps(2);
    tdeps[0] = GetRawPtr(tangential_x);
    tdeps[1] = GetRawPtr(W);
    std::vector<Number> sdeps(1);
    sdeps[0] = pd_pert_x;

    if (!curr_Wu_x_cache_.GetCachedResult(result, tdeps, sdeps)) {
      result = curr_W_times_vec_x(*tangential_x);
      curr_Wu_x_cache_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;
  }

  SmartPtr<const Vector>
  InexactCq::curr_Wu_s()
  {
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> tangential_s = InexData().tangential_s();
    SmartPtr<const Vector> sigma_s = ip_cq_->curr_sigma_s();
    Number pd_pert_x;
    Number pd_pert_s;
    Number pd_pert_c;
    Number pd_pert_d;
    ip_data_->getPDPert(pd_pert_x, pd_pert_s, pd_pert_c, pd_pert_d);

    std::vector<const TaggedObject*> tdeps(2);
    tdeps[0] = GetRawPtr(tangential_s);
    tdeps[1] = GetRawPtr(sigma_s);
    std::vector<Number> sdeps(1);
    sdeps[0] = pd_pert_s;

    if (!curr_Wu_s_cache_.GetCachedResult(result, tdeps, sdeps)) {
      result = curr_W_times_vec_s(*tangential_s);
      curr_Wu_s_cache_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;
  }

  Number
  InexactCq::curr_uWu()
  {
    Number result;

    SmartPtr<const Vector> tangential_x = InexData().tangential_x();
    SmartPtr<const Vector> tangential_s = InexData().tangential_s();
    SmartPtr<const Vector> Wu_x = curr_Wu_x();
    SmartPtr<const Vector> Wu_s = curr_Wu_s();

    std::vector<const TaggedObject*> tdeps(4);
    tdeps[0] = GetRawPtr(tangential_x);
    tdeps[1] = GetRawPtr(tangential_s);
    tdeps[2] = GetRawPtr(Wu_x);
    tdeps[3] = GetRawPtr(Wu_s);

    if (!curr_uWu_cache_.GetCachedResult(result, tdeps)) {
      result = Wu_x->Dot(*tangential_x) + Wu_s->Dot(*tangential_s);
      curr_uWu_cache_.AddCachedResult(result, tdeps);
    }

    return result;
  }

  SmartPtr<const Vector>
  InexactCq::curr_jac_times_normal_c()
  {
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> normal_x = InexData().normal_x();

    if (!curr_jac_times_normal_c_cache_.GetCachedResult2Dep(result, *x, *normal_x)) {
      result = ip_cq_->curr_jac_c_times_vec(*normal_x);
      curr_jac_times_normal_c_cache_.AddCachedResult2Dep(result, *x, *normal_x);
    }

    return result;
  }

  SmartPtr<const Vector>
  InexactCq::curr_jac_times_normal_d()
  {
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> normal_x = InexData().normal_x();
    SmartPtr<const Vector> normal_s = InexData().normal_s();

    if (!curr_jac_times_normal_d_cache_.GetCachedResult3Dep(result, *x, *normal_x, *normal_s)) {
      SmartPtr<const Vector> Ax = ip_cq_->curr_jac_d_times_vec(*normal_x);
      SmartPtr<Vector> tmp = Ax->MakeNew();
      tmp->AddTwoVectors(1., *Ax, -1., *normal_s, 0.);

      result = ConstPtr(tmp);
      curr_jac_times_normal_d_cache_.AddCachedResult3Dep(result, *x, *normal_x, *normal_s);
    }

    return result;
  }

} // namespace Ipopt
