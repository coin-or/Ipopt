// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter              IBM    2004-09-23

#include "IpWarmStartIterateInitializer.hpp"

// DELETEME
#include "IpDenseVector.hpp"
#ifdef OLD_C_HEADERS
# include <math.h>
#else
# include <cmath>
#endif

namespace Ipopt
{
  DBG_SET_VERBOSITY(0);

  WarmStartIterateInitializer::WarmStartIterateInitializer()
      :
      IterateInitializer()
  {}

  bool WarmStartIterateInitializer::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    Number value = 0.0;
    // Check for the algorithm options

    if (options.GetNumericValue("warm_start_bound_push", value, prefix)) {
      ASSERT_EXCEPTION(value > 0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"warm_start_bound_push\": This value must be larger than 0.");
      warm_start_bound_push_ = value;
    }
    else {
      warm_start_bound_push_ = 1e-6;
    }

    if (options.GetNumericValue("warm_start_bound_frac", value, prefix)) {
      ASSERT_EXCEPTION(value > 0 && value <= 0.5, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"warm_start_bound_frac\": Value must be between 0 and 0.5.");
      warm_start_bound_frac_ = value;
    }
    else {
      warm_start_bound_frac_ = 1e-6;
    }

    if (options.GetNumericValue("warm_start_mult_bound_push", value, prefix)) {
      ASSERT_EXCEPTION(value > 0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"warm_start_mult_bound_push\": This value must be larger than 0.");
      warm_start_mult_bound_push_ = value;
    }
    else {
      warm_start_mult_bound_push_ = 1e-6;
    }

    if (options.GetNumericValue("warm_start_mult_init_max", value, prefix)) {
      warm_start_mult_init_max_ = value;
    }
    else {
      warm_start_mult_init_max_ = 1e6;
    }

    if (options.GetNumericValue("warm_start_target_mu", value, prefix)) {
      warm_start_target_mu_ = value;
    }
    else {
      warm_start_target_mu_ = 1e-3;
    }

    return true;
  }

  bool WarmStartIterateInitializer::SetInitialIterates()
  {
    DBG_START_METH("WarmStartIterateInitializer::SetInitialIterates",
                   dbg_verbosity);

    // Get the starting values provided by the NLP and store them
    //.in the ip_data current fields.  The following line only requests
    // intial values for the primal variables x, but later we might
    // make this more flexible based on user options.

    /////////////////////////////////////////////////////////////////////
    //                   Initialize primal variables                   //
    /////////////////////////////////////////////////////////////////////

    // Get the intial values for x, y_c, y_d, z_L, z_U,
    IpData().InitializeDataStructures(IpNLP(), true, true, true,
                                      true, true, false, false);

    Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION, "user-provided x",
                        *IpData().curr_x());
    Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION, "user-provided y_c",
                        *IpData().curr_y_c());
    Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION, "user-provided y_d",
                        *IpData().curr_y_d());
    Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION, "user-provided z_L",
                        *IpData().curr_z_L());
    Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION, "user-provided z_U",
                        *IpData().curr_z_U());
    if (Jnlst().ProduceOutput(J_MOREVECTOR, J_INITIALIZATION)) {
      Jnlst().PrintVector(J_MOREVECTOR, J_INITIALIZATION, "d at user-provided x",
                          *IpCq().curr_d());
    }

    SmartPtr<Vector> tmp;

    {
      // If requested, make sure that the multipliers are not too large
      SmartPtr<const Vector> new_y_c;
      SmartPtr<const Vector> new_y_d;
      SmartPtr<const Vector> new_z_L;
      SmartPtr<const Vector> new_z_U;
      if (warm_start_mult_init_max_>0.) {
        SmartPtr<Vector> y_c = IpData().curr_y_c()->MakeNew();
        y_c->Copy(*IpData().curr_y_c());
        tmp = y_c->MakeNew();
        tmp->Set(warm_start_mult_init_max_);
        y_c->ElementWiseMin(*tmp);
        tmp->Set(-warm_start_mult_init_max_);
        y_c->ElementWiseMax(*tmp);
        new_y_c = ConstPtr(y_c);

        SmartPtr<Vector> y_d = IpData().curr_y_d()->MakeNew();
        y_d->Copy(*IpData().curr_y_d());
        tmp = y_d->MakeNew();
        tmp->Set(warm_start_mult_init_max_);
        y_d->ElementWiseMin(*tmp);
        tmp->Set(-warm_start_mult_init_max_);
        y_d->ElementWiseMax(*tmp);
        new_y_d = ConstPtr(y_d);

        SmartPtr<Vector> z_L = IpData().curr_z_L()->MakeNew();
        z_L->Copy(*IpData().curr_z_L());
        tmp = z_L->MakeNew();
        tmp->Set(warm_start_mult_init_max_);
        z_L->ElementWiseMin(*tmp);
        new_z_L = ConstPtr(z_L);

        SmartPtr<Vector> z_U = IpData().curr_z_U()->MakeNew();
        z_U->Copy(*IpData().curr_z_U());
        tmp = z_U->MakeNew();
        tmp->Set(warm_start_mult_init_max_);
        z_U->ElementWiseMin(*tmp);
        new_z_U = ConstPtr(z_U);
      }
      else {
        new_y_c = IpData().curr_y_c();
        new_y_d = IpData().curr_y_d();
        new_z_L = IpData().curr_z_L();
        new_z_U = IpData().curr_z_U();
      }

      // Get the initial values for v_L and v_U out of y_d
      SmartPtr<Vector> new_v_L = IpNLP().d_L()->MakeNew();
      IpNLP().Pd_L()->TransMultVector(-1., *new_y_d, 0., *new_v_L);
      tmp = new_v_L->MakeNew();
      tmp->Set(warm_start_mult_bound_push_);
      new_v_L->ElementWiseMax(*tmp);

      SmartPtr<Vector> new_v_U = IpNLP().d_U()->MakeNew();
      IpNLP().Pd_U()->TransMultVector(1., *new_y_d, 0., *new_v_U);
      tmp = new_v_U->MakeNew();
      tmp->Set(warm_start_mult_bound_push_);
      new_v_U->ElementWiseMax(*tmp);

      // Make those currected values current (and initialize s)
      IpData().SetTrialPrimalVariablesFromPtr(IpData().curr_x(),
                                              IpCq().curr_d());
      IpData().SetTrialConstraintMultipliersFromPtr(new_y_c, new_y_d);
      IpData().SetTrialBoundMultipliersFromPtr(new_z_L, new_z_U,
          ConstPtr(new_v_L),
          ConstPtr(new_v_U));
      IpData().AcceptTrialPoint();
    }

    // Now apply the target mu heuristic if required
    if (warm_start_target_mu_>0.) {
      SmartPtr<Vector> new_x;
      SmartPtr<Vector> new_z_L;
      process_target_mu(1., *IpData().curr_x(), *IpCq().curr_slack_x_L(),
                        *IpData().curr_z_L(), *IpNLP().Px_L(),
                        new_x, new_z_L);
      SmartPtr<Vector> new_s;
      SmartPtr<Vector> new_v_L;
      process_target_mu(1., *IpData().curr_s(), *IpCq().curr_slack_s_L(),
                        *IpData().curr_v_L(), *IpNLP().Pd_L(),
                        new_s, new_v_L);
      IpData().SetTrialPrimalVariablesFromPtr(ConstPtr(new_x),
                                              ConstPtr(new_s));

      SmartPtr<Vector> new_z_U;
      process_target_mu(-1., *IpData().trial_x(), *IpCq().trial_slack_x_U(),
                        *IpData().curr_z_U(), *IpNLP().Px_U(),
                        new_x, new_z_U);
      SmartPtr<Vector> new_v_U;
      process_target_mu(-1., *IpData().trial_s(), *IpCq().trial_slack_s_U(),
                        *IpData().curr_v_U(), *IpNLP().Pd_U(),
                        new_s, new_v_U);
      Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION, "DELETEME new_s",
                          *new_s);
      Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION, "DELETEME d_U",
                          *IpNLP().d_U());
      Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION, "DELETEME new_v_U",
                          *new_v_U);


      IpData().SetTrialPrimalVariablesFromPtr(ConstPtr(new_x),
                                              ConstPtr(new_s));
      IpData().SetTrialConstraintMultipliersFromPtr(IpData().curr_y_c(),
          IpData().curr_y_d());
      IpData().SetTrialBoundMultipliersFromPtr(ConstPtr(new_z_L),
          ConstPtr(new_z_U),
          ConstPtr(new_v_L),
          ConstPtr(new_v_U));
      IpData().AcceptTrialPoint();

      // We need to call this to make sure that we don't get an error
      // message because at some point a slack became too small
      IpCq().ResetAdjustedTrialSlacks();
    }
    Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION, "DELETEME curr_s",
                        *IpData().curr_s());

    {
      SmartPtr<const Vector> new_x;
      SmartPtr<const Vector> new_s;
      // Push the primal x variables
      push_variables(warm_start_bound_push_,
                     warm_start_bound_frac_,
                     "x",
                     *IpData().curr_x(),
                     new_x,
                     *IpNLP().x_L(),
                     *IpNLP().x_U(),
                     *IpNLP().Px_L(),
                     *IpNLP().Px_U());

      IpData().SetTrialPrimalVariablesFromPtr(new_x, new_s);

      // Push the primal s variables
      push_variables(warm_start_bound_push_,
                     warm_start_bound_frac_,
                     "s",
                     *IpData().curr_s(),
                     new_s,
                     *IpNLP().d_L(),
                     *IpNLP().d_U(),
                     *IpNLP().Pd_L(),
                     *IpNLP().Pd_U());
      Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION, "DELETEME new_s cor",
                          *new_s);

      // Push the multipliers
      SmartPtr<Vector> new_z_L = IpData().curr_z_L()->MakeNew();
      new_z_L->Copy(*IpData().curr_z_L());
      tmp = IpData().curr_z_L()->MakeNew();
      tmp->Set(warm_start_mult_bound_push_);
      new_z_L->ElementWiseMax(*tmp);

      SmartPtr<Vector> new_z_U = IpData().curr_z_U()->MakeNew();
      new_z_U->Copy(*IpData().curr_z_U());
      tmp = IpData().curr_z_U()->MakeNew();
      tmp->Set(warm_start_mult_bound_push_);
      new_z_U->ElementWiseMax(*tmp);

      SmartPtr<Vector> new_v_L = IpData().curr_v_L()->MakeNew();
      new_v_L->Copy(*IpData().curr_v_L());
      tmp = IpData().curr_v_L()->MakeNew();
      tmp->Set(warm_start_mult_bound_push_);
      new_v_L->ElementWiseMax(*tmp);

      SmartPtr<Vector> new_v_U = IpData().curr_v_U()->MakeNew();
      new_v_U->Copy(*IpData().curr_v_U());
      tmp = IpData().curr_v_U()->MakeNew();
      tmp->Set(warm_start_mult_bound_push_);
      new_v_U->ElementWiseMax(*tmp);

      // Make sure the new variables are current
      IpData().SetTrialPrimalVariablesFromPtr(new_x, new_s);
      IpData().SetTrialConstraintMultipliersFromPtr(IpData().curr_y_c(),
          IpData().curr_y_d());
      IpData().SetTrialBoundMultipliersFromPtr(ConstPtr(new_z_L),
          ConstPtr(new_z_U),
          ConstPtr(new_v_L),
          ConstPtr(new_v_U));
      IpData().AcceptTrialPoint();
    }

    Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION, "initial x",
                        *IpData().curr_x());
    Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION, "initial s",
                        *IpData().curr_s());
    Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION, "initial y_c",
                        *IpData().curr_y_c());
    Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION, "initial y_d",
                        *IpData().curr_y_d());
    Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION, "initial z_L",
                        *IpData().curr_z_L());
    Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION, "initial z_U",
                        *IpData().curr_z_U());
    Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION, "initial v_L",
                        *IpData().curr_v_L());
    Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION, "initial v_U",
                        *IpData().curr_v_U());
    if (Jnlst().ProduceOutput(J_MOREVECTOR, J_INITIALIZATION)) {
      Jnlst().PrintVector(J_MOREVECTOR, J_INITIALIZATION, "initial slack_x_L",
                          *IpCq().curr_slack_x_L());
      Jnlst().PrintVector(J_MOREVECTOR, J_INITIALIZATION, "initial slack_x_U",
                          *IpCq().curr_slack_x_U());
      Jnlst().PrintVector(J_MOREVECTOR, J_INITIALIZATION, "initial slack_s_L",
                          *IpCq().curr_slack_s_L());
      Jnlst().PrintVector(J_MOREVECTOR, J_INITIALIZATION, "initial slack_s_U",
                          *IpCq().curr_slack_s_U());
    }

    return true;
  }

  void WarmStartIterateInitializer::process_target_mu(Number factor,
      const Vector& curr_vars,
      const Vector& curr_slacks,
      const Vector& curr_mults,
      const Matrix& P,
      SmartPtr<Vector>& new_vars,
      SmartPtr<Vector>& new_mults)
  {
    SmartPtr<Vector> new_slacks = curr_slacks.MakeNew();
    new_slacks->Copy(curr_slacks);
    new_mults = curr_mults.MakeNew();
    new_mults->Copy(curr_mults);
    adapt_to_target_mu(*new_slacks, *new_mults, warm_start_target_mu_);
    Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION, "DELETEME new_slacks",
                        *new_slacks);
    new_slacks->Axpy(-1, curr_slacks); // this is now correction step
    new_vars = curr_vars.MakeNew();
    new_vars->Copy(curr_vars);
    P.MultVector(factor, *new_slacks, 1., *new_vars);
  }

  void WarmStartIterateInitializer::push_variables(Number bound_bush,
      Number bound_frac,
      std::string name,
      const Vector& orig_x,
      SmartPtr<const Vector>& new_x,
      const Vector& x_L,
      const Vector& x_U,
      const Matrix& Px_L,
      const Matrix& Px_U)
  {
    DBG_START_METH("WarmStartIterateInitializer::push_variables",
                   dbg_verbosity);
    // Calculate any required shift in x0 and s0
    const double dbl_min = std::numeric_limits<double>::min();
    const double tiny_double = 100.0*dbl_min;

    SmartPtr<Vector> tmp = orig_x.MakeNew();
    SmartPtr<Vector> tmp_l = x_L.MakeNew();
    SmartPtr<Vector> tmp_u = x_U.MakeNew();
    SmartPtr<Vector> tiny_l = x_L.MakeNew();
    tiny_l->Set(tiny_double);

    // Calculate p_l
    SmartPtr<Vector> q_l = x_L.MakeNew();
    SmartPtr<Vector> p_l = x_L.MakeNew();

    DBG_PRINT_MATRIX(2,"Px_L", Px_L);
    DBG_PRINT_VECTOR(2, "x_L", x_L);
    DBG_PRINT_VECTOR(2, "tmp", *tmp);
    Px_L.MultVector(1.0, x_L, 0.0, *tmp);
    Px_U.TransMultVector(1.0, *tmp, 0.0, *tmp_u);
    tmp_u->AddOneVector(1., x_U, -1.);
    Px_U.MultVector(1.0, *tmp_u, 0.0, *tmp);
    Px_L.TransMultVector(1.0, *tmp, 0.0, *q_l);
    q_l->AddOneVector(-1.0, *tiny_l, warm_start_bound_frac_);

    tmp_l->Set(1.0);
    p_l->Copy(x_L);
    p_l->ElementWiseSgn();
    p_l->ElementWiseMultiply(x_L);
    p_l->ElementWiseMax(*tmp_l);
    p_l->AddOneVector(-1.0, *tiny_l, warm_start_bound_push_);

    q_l->ElementWiseReciprocal();
    p_l->ElementWiseReciprocal();

    p_l->ElementWiseMax(*q_l);
    p_l->ElementWiseReciprocal();
    p_l->Axpy(1.0, *tiny_l);

    // Calculate p_u
    SmartPtr<Vector> q_u = x_U.MakeNew();
    SmartPtr<Vector> p_u = x_U.MakeNew();
    SmartPtr<Vector> tiny_u = x_U.MakeNew();
    tiny_u->Set(tiny_double);

    Px_U.MultVector(1.0, x_U, 0.0, *tmp);
    Px_L.TransMultVector(1.0, *tmp, 0.0, *tmp_l);
    tmp_l->Axpy(-1.0, x_L);
    Px_L.MultVector(1.0, *tmp_l, 0.0, *tmp);
    Px_U.TransMultVector(1.0, *tmp, 0.0, *q_u);
    q_u->AddOneVector(-1.0, *tiny_u, warm_start_bound_frac_);
    DBG_PRINT_VECTOR(2,"q_u",*q_u);

    tmp_u->Set(1.0);
    p_u->Copy(x_U);
    p_u->ElementWiseSgn();
    p_u->ElementWiseMultiply(x_U);
    p_u->ElementWiseMax(*tmp_u);
    p_u->AddOneVector(-1.0, *tiny_u, warm_start_bound_push_);
    DBG_PRINT_VECTOR(2,"p_u",*p_u);

    q_u->ElementWiseReciprocal();
    p_u->ElementWiseReciprocal();

    p_u->ElementWiseMax(*q_u);
    p_u->ElementWiseReciprocal();
    p_u->Axpy(1.0, *tiny_u);
    DBG_PRINT_VECTOR(2,"actual_p_u",*p_u);

    // Calculate the new x
    SmartPtr<Vector> delta_x = orig_x.MakeNew();

    SmartPtr<Vector> zero_l = x_L.MakeNew();
    zero_l->Set(0.0);
    SmartPtr<Vector> zero_u = x_U.MakeNew();
    zero_u->Set(0.0);

    Px_L.TransMultVector(-1.0, orig_x, 0.0, *tmp_l);
    tmp_l->AddTwoVectors(1.0, x_L, 1.0, *p_l, 1.);
    tmp_l->ElementWiseMax(*zero_l);
    Number nrm_l = tmp_l->Amax();
    if (nrm_l>0.) {
      Px_L.MultVector(1.0, *tmp_l, 0.0, *delta_x);
    }
    else {
      delta_x->Set(0.);
    }

    Px_U.TransMultVector(1.0, orig_x, 0.0, *tmp_u);
    tmp_u->AddTwoVectors(-1.0, x_U, 1.0, *p_u, 1.);
    tmp_u->ElementWiseMax(*zero_u);
    Number nrm_u = tmp_u->Amax();
    if (nrm_u>0.) {
      Px_U.MultVector(-1.0, *tmp_u, 1.0, *delta_x);
    }

    if (nrm_l > 0 || nrm_u > 0) {
      delta_x->Axpy(1.0, orig_x);
      new_x = ConstPtr(delta_x);
      Jnlst().Printf(J_DETAILED, J_INITIALIZATION, "Moved initial values of %s sufficiently inside the bounds.\n", name.c_str());
      Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION, "original vars", orig_x);
      Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION, "new vars", *new_x);
    }
    else {
      new_x = &orig_x;
    }
  }

  void WarmStartIterateInitializer::adapt_to_target_mu(Vector& new_s,
      Vector& new_z,
      Number target_mu)
  {
    DBG_ASSERT(new_s.Dim() == new_z.Dim());

    DenseVector* dnew_s = dynamic_cast<DenseVector*>(&new_s);
    assert(dnew_s);
    DenseVector* dnew_z = dynamic_cast<DenseVector*>(&new_z);
    assert(dnew_z);
    Number* values_s = dnew_s->Values();
    Number* values_z = dnew_z->Values();

    for (Index i=0; i<new_s.Dim(); i++) {
      if (values_s[i] > 1e4*values_z[i]) {
        values_z[i] = target_mu/values_s[i];
        if (values_z[i]>values_s[i]) {
          values_s[i] = values_z[i] = sqrt(target_mu);
        }
      }
      else if (values_z[i] > 1e4*values_s[i]) {
        values_s[i] = target_mu/values_z[i];
        if (values_s[i]>values_z[i]) {
          values_s[i] = values_z[i] = sqrt(target_mu);
        }
      }
      else {
        values_s[i] = values_z[i] = sqrt(target_mu);
      }
    }
  }

} // namespace Ipopt
