// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter              IBM    2004-09-23

#include "IpDefaultIterateInitializer.hpp"

namespace Ipopt
{
  DBG_SET_VERBOSITY(0);

  DefaultIterateInitializer::DefaultIterateInitializer
  (const SmartPtr<EqMultiplierCalculator>& eq_mult_calculator)
      :
      IterateInitializer(),
      eq_mult_calculator_(eq_mult_calculator)
  {}

  bool DefaultIterateInitializer::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    Number value = 0.0;
    // Check for the algorithm options
    if (options.GetNumericValue("bound_push", value, prefix)) {
      ASSERT_EXCEPTION(value > 0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"bound_push\": This value must be larger than 0.");
      bound_push_ = value;
    }
    else {
      bound_push_ = 0.01;
    }

    if (options.GetNumericValue("bound_frac", value, prefix)) {
      ASSERT_EXCEPTION(value > 0 && value < 0.5, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"bound_frac\": Value must be between 0 and 0.5.");
      bound_frac_ = value;
    }
    else {
      bound_frac_ = 0.01;
    }

    if (options.GetNumericValue("laminitmax", value, prefix)) {
      ASSERT_EXCEPTION(value >= 0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"laminitmax\": Value must be non-negative.");
      laminitmax_ = value;
    }
    else {
      laminitmax_ = 1e3;
    }

    if (options.GetNumericValue("boundmultinitval", value, prefix)) {
      ASSERT_EXCEPTION(value > 0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"boundmultinitval\": Value must be positive.");
      boundmultinitval_ = value;
    }
    else {
      boundmultinitval_ = 1.;
    }

    bool retvalue = true;
    if (IsValid(eq_mult_calculator_)) {
      retvalue = eq_mult_calculator_->Initialize(Jnlst(), IpNLP(), IpData(),
                 IpCq(), options, prefix);
    }
    return retvalue;
  }

  bool DefaultIterateInitializer::SetInitialIterates()
  {
    DBG_START_METH("DefaultIterateInitializer::SetInitialIterates",
                   dbg_verbosity);

    // Get the starting values provided by the NLP and store them
    //.in the ip_data current fields.  The following line only requests
    // intial values for the primal variables x, but later we might
    // make this more flexible based on user options.

    /////////////////////////////////////////////////////////////////////
    //                   Initialize primal variables                   //
    /////////////////////////////////////////////////////////////////////

    IpData().InitializeDataStructures(IpNLP(), true, false, false,
                                      false, false, false, false);

    DBG_PRINT_VECTOR(2, "curr_x", *IpData().curr_x());

    // Now we compute the initial values that the algorithm is going to
    // actually use.  We first store them in the trial fields in ip_data.

    // Calculate any required shift in x0 and s0
    const double dbl_min = std::numeric_limits<double>::min();
    const double tiny_double = 100.0*dbl_min;
    SmartPtr<const Vector> x = IpData().curr_x();
    SmartPtr<const Vector> x_L = IpNLP().x_L();
    SmartPtr<const Vector> x_U = IpNLP().x_U();
    SmartPtr<const Matrix> Px_L = IpNLP().Px_L();
    SmartPtr<const Matrix> Px_U = IpNLP().Px_U();

    SmartPtr<Vector> tmp = x->MakeNew();
    SmartPtr<Vector> tmp_l = x_L->MakeNew();
    SmartPtr<Vector> tmp_u = x_U->MakeNew();
    SmartPtr<Vector> tiny_l = x_L->MakeNew();
    tiny_l->Set(tiny_double);

    // Calculate p_l
    SmartPtr<Vector> q_l = x_L->MakeNew();
    SmartPtr<Vector> p_l = x_L->MakeNew();

    DBG_PRINT_MATRIX(2,"Px_L", *Px_L);
    DBG_PRINT_VECTOR(2, "x_L", *x_L);
    DBG_PRINT_VECTOR(2, "tmp", *tmp);
    Px_L->MultVector(1.0, *x_L, 0.0, *tmp);
    Px_U->TransMultVector(1.0, *tmp, 0.0, *tmp_u);
    tmp_u->AddOneVector(1., *x_U, -1.);
    /* DELE
    tmp_u->Axpy(-1.0, *x_U);
    tmp_u->Scal(-1.0);
    */
    Px_U->MultVector(1.0, *tmp_u, 0.0, *tmp);
    Px_L->TransMultVector(1.0, *tmp, 0.0, *q_l);
    q_l->AddOneVector(-1.0, *tiny_l, bound_frac_);
    /* DELE
    q_l->Scal(bound_frac_);
    q_l->Axpy(-1.0, *tiny_l);
    */

    tmp_l->Set(1.0);
    p_l->Copy(*x_L);
    p_l->ElementWiseSgn();
    p_l->ElementWiseMultiply(*x_L);
    p_l->ElementWiseMax(*tmp_l);
    p_l->AddOneVector(-1.0, *tiny_l, bound_push_);
    /* DELE
    p_l->Scal(bound_push_);
    p_l->Axpy(-1.0, *tiny_l);
    */

    q_l->ElementWiseReciprocal();
    p_l->ElementWiseReciprocal();

    p_l->ElementWiseMax(*q_l);
    p_l->ElementWiseReciprocal();
    p_l->Axpy(1.0, *tiny_l);

    // Calculate p_u
    SmartPtr<Vector> q_u = x_U->MakeNew();
    SmartPtr<Vector> p_u = x_U->MakeNew();
    SmartPtr<Vector> tiny_u = x_U->MakeNew();
    tiny_u->Set(tiny_double);

    Px_U->MultVector(1.0, *x_U, 0.0, *tmp);
    Px_L->TransMultVector(1.0, *tmp, 0.0, *tmp_l);
    tmp_l->Axpy(-1.0, *x_L);
    Px_L->MultVector(1.0, *tmp_l, 0.0, *tmp);
    Px_U->TransMultVector(1.0, *tmp, 0.0, *q_u);
    q_u->AddOneVector(-1.0, *tiny_u, bound_frac_);
    /* DELE
    q_u->Scal(bound_frac_);
    q_u->Axpy(-1.0, *tiny_u);
    */
    DBG_PRINT_VECTOR(2,"q_u",*q_u);

    tmp_u->Set(1.0);
    p_u->Copy(*x_U);
    p_u->ElementWiseSgn();
    p_u->ElementWiseMultiply(*x_U);
    p_u->ElementWiseMax(*tmp_u);
    p_u->AddOneVector(-1.0, *tiny_u, bound_push_);
    /* DELE
    p_u->Scal(bound_push_);
    p_u->Axpy(-1.0, *tiny_u);
    */
    DBG_PRINT_VECTOR(2,"p_u",*p_u);


    q_u->ElementWiseReciprocal();
    p_u->ElementWiseReciprocal();

    p_u->ElementWiseMax(*q_u);
    p_u->ElementWiseReciprocal();
    p_u->Axpy(1.0, *tiny_u);
    DBG_PRINT_VECTOR(2,"actual_p_u",*p_u);


    // Calculate the new x
    SmartPtr<Vector> delta_x = x->MakeNew();

    SmartPtr<Vector> zero_l = x_L->MakeNew();
    zero_l->Set(0.0);
    SmartPtr<Vector> zero_u = x_U->MakeNew();
    zero_u->Set(0.0);

    Px_L->TransMultVector(-1.0, *x, 0.0, *tmp_l);
    tmp_l->AddTwoVectors(1.0, *x_L, 1.0, *p_l, 1.);
    /* DELE
    tmp_l->Axpy(1.0, *x_L);
    tmp_l->Axpy(1.0, *p_l);
    */
    tmp_l->ElementWiseMax(*zero_l);
    Number nrm_l = tmp_l->Amax();
    if (nrm_l>0.) {
      Px_L->MultVector(1.0, *tmp_l, 0.0, *delta_x);
    }
    else {
      delta_x->Set(0.);
    }

    Px_U->TransMultVector(1.0, *x, 0.0, *tmp_u);
    tmp_u->AddTwoVectors(-1.0, *x_U, 1.0, *p_u, 1.);
    /* DELE
    tmp_u->Axpy(-1.0, *x_U);
    tmp_u->Axpy(1.0, *p_u);
    */
    tmp_u->ElementWiseMax(*zero_u);
    Number nrm_u = tmp_u->Amax();
    if (nrm_u>0.) {
      Px_U->MultVector(-1.0, *tmp_u, 1.0, *delta_x);
    }
    //delta_x->Axpy(-1.0, *tmp);

    SmartPtr<Vector> new_x = delta_x;
    new_x->Axpy(1.0, *x);

    IpData().SetTrialXVariables(*new_x);
    if (nrm_l > 0 || nrm_u > 0) {
      Jnlst().Printf(J_DETAILED, J_INITIALIZATION, "Moved initial values of x sufficiently inside the bounds.\n");
      Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION, "original x", *IpData().curr_x());
      Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION, "new x", *IpData().trial_x());
    }

    // Calculate the shift in s...
    SmartPtr<const Vector> s = IpCq().trial_d();
    DBG_PRINT_VECTOR(2, "s", *s);

    SmartPtr<const Vector> d_L = IpNLP().d_L();
    SmartPtr<const Vector> d_U = IpNLP().d_U();
    SmartPtr<const Matrix> Pd_L = IpNLP().Pd_L();
    SmartPtr<const Matrix> Pd_U = IpNLP().Pd_U();
    DBG_PRINT_VECTOR(2, "d_L", *d_L);
    DBG_PRINT_VECTOR(2, "d_U", *d_U);

    // bump the starting point if necessary
    tmp = s->MakeNew();
    tmp_l = d_L->MakeNew();
    tmp_u = d_U->MakeNew();
    tiny_l = d_L->MakeNew();
    tiny_l->Set(tiny_double);

    // Calculate p_l
    q_l = d_L->MakeNew();
    p_l = d_L->MakeNew();

    Pd_L->MultVector(1.0, *d_L, 0.0, *tmp);
    Pd_U->TransMultVector(1.0, *tmp, 0.0, *tmp_u);
    tmp_u->AddOneVector(1., *d_U, -1.);
    /* DELE
    tmp_u->Axpy(-1.0, *d_U);
    tmp_u->Scal(-1.0);
    */
    Pd_U->MultVector(1.0, *tmp_u, 0.0, *tmp);
    Pd_L->TransMultVector(1.0, *tmp, 0.0, *q_l);
    q_l->AddOneVector(-1.0, *tiny_l, bound_frac_);
    /* DELE
    q_l->Scal(bound_frac_);
    q_l->Axpy(-1.0, *tiny_l);
    */

    tmp_l->Set(1.0);
    p_l->Copy(*d_L);
    p_l->ElementWiseSgn();
    p_l->ElementWiseMultiply(*d_L);
    p_l->ElementWiseMax(*tmp_l);
    p_l->AddOneVector(-1.0, *tiny_l, bound_push_);
    /* DELE
    p_l->Scal(bound_push_);
    p_l->Axpy(-1.0, *tiny_l);
    */

    q_l->ElementWiseReciprocal();
    p_l->ElementWiseReciprocal();

    p_l->ElementWiseMax(*q_l);
    p_l->ElementWiseReciprocal();
    p_l->Axpy(1.0, *tiny_l);
    DBG_PRINT_VECTOR(2, "p_l", *p_l);

    // Calculate p_u
    q_u = d_U->MakeNew();
    p_u = d_U->MakeNew();
    tiny_u = d_U->MakeNew();
    tiny_u->Set(tiny_double);

    Pd_U->MultVector(1.0, *d_U, 0.0, *tmp);
    Pd_L->TransMultVector(1.0, *tmp, 0.0, *tmp_l);
    tmp_l->Axpy(-1.0, *d_L);
    Pd_L->MultVector(1.0, *tmp_l, 0.0, *tmp);
    Pd_U->TransMultVector(1.0, *tmp, 0.0, *q_u);
    q_u->AddOneVector(-1.0, *tiny_u, bound_frac_);
    /* DELE
    q_u->Scal(bound_frac_);
    q_u->Axpy(-1.0, *tiny_u);
    */

    tmp_u->Set(1.0);
    p_u->Copy(*d_U);
    p_u->ElementWiseSgn();
    p_u->ElementWiseMultiply(*d_U);
    p_u->ElementWiseMax(*tmp_u);
    p_u->AddOneVector(-1.0, *tiny_u, bound_push_);
    /* DELE
    p_u->Scal(bound_push_);
    p_u->Axpy(-1.0, *tiny_u);
    */

    q_u->ElementWiseReciprocal();
    p_u->ElementWiseReciprocal();

    p_u->ElementWiseMax(*q_u);
    p_u->ElementWiseReciprocal();
    p_u->Axpy(1.0, *tiny_u);
    DBG_PRINT_VECTOR(2, "p_u", *p_u);

    zero_l = d_L->MakeNew();
    zero_l->Set(0.0);
    zero_u = d_U->MakeNew();
    zero_u->Set(0.0);

    Pd_L->TransMultVector(-1.0, *s, 0.0, *tmp_l);
    tmp_l->AddTwoVectors(1.0, *d_L, 1.0, *p_l, 1.);
    /* DELE
    tmp_l->Axpy(1.0, *d_L);
    tmp_l->Axpy(1.0, *p_l);
    */
    tmp_l->ElementWiseMax(*zero_l);
    nrm_l = tmp_l->Amax();
    SmartPtr<Vector> delta_s = s->MakeNew();
    if (nrm_l>0.) {
      Pd_L->MultVector(1.0, *tmp_l, 0.0, *delta_s);
    }
    else {
      delta_s->Set(0.);
    }

    Pd_U->TransMultVector(1.0, *s, 0.0, *tmp_u);
    tmp_u->AddTwoVectors(-1.0, *d_U, 1.0, *p_u, 1.);
    /* DELE
    tmp_u->Axpy(-1.0, *d_U);
    tmp_u->Axpy(1.0, *p_u);
    */
    tmp_u->ElementWiseMax(*zero_u);
    nrm_u = tmp_u->Amax();
    Pd_U->MultVector(-1.0, *tmp_u, 1.0, *delta_s);

    SmartPtr<Vector> new_s = delta_s;
    new_s->Axpy(1.0, *s);

    IpData().SetTrialSVariables(*new_s);
    if (nrm_l > 0 || nrm_u > 0) {
      Jnlst().Printf(J_DETAILED, J_INITIALIZATION,
                     "Moved initial values of s sufficiently inside the bounds.\n");
      Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION,
                          "original s", *s);
      Jnlst().PrintVector(J_VECTOR, J_INITIALIZATION,
                          "new s", *IpData().trial_s());
    }

    /////////////////////////////////////////////////////////////////////
    //                   Initialize bound multipliers                  //
    /////////////////////////////////////////////////////////////////////

    // Initialize the bound multipliers to 1.
    SmartPtr<Vector> z_L = IpData().curr_z_L()->MakeNew();
    SmartPtr<Vector> z_U = IpData().curr_z_U()->MakeNew();
    SmartPtr<Vector> v_L = IpData().curr_v_L()->MakeNew();
    SmartPtr<Vector> v_U = IpData().curr_v_U()->MakeNew();
    z_L->Set(boundmultinitval_);  //TODO make this a parameter
    z_U->Set(boundmultinitval_);
    v_L->Set(boundmultinitval_);
    v_U->Set(boundmultinitval_);
    IpData().SetTrialBoundMultipliers(*z_L, *z_U, *v_L, *v_U);

    /////////////////////////////////////////////////////////////////////
    //           Initialize equality constraint multipliers            //
    /////////////////////////////////////////////////////////////////////

    if (IsValid(eq_mult_calculator_) && laminitmax_>0.) {
      // First move all the trial data into the current fields, since
      // those values are needed to compute the initial values for
      // the multipliers
      IpData().CopyTrialToCurrent();
      SmartPtr<Vector> y_c = IpData().curr_y_c()->MakeNew();
      SmartPtr<Vector> y_d = IpData().curr_y_d()->MakeNew();
      bool retval = eq_mult_calculator_->CalculateMultipliers(*y_c, *y_d);
      if (!retval) {
        y_c->Set(0.0);
        y_d->Set(0.0);
      }
      else {
        Jnlst().Printf(J_DETAILED, J_INITIALIZATION,
                       "Least square estimates max(y_c) = %e, max(y_d) = %e\n",
                       y_c->Amax(), y_d->Amax());
        Number laminitnrm = Max(y_c->Amax(), y_d->Amax());
        if (laminitnrm > laminitmax_) {
          y_c->Set(0.0);
          y_d->Set(0.0);
        }
      }
      IpData().SetTrialEqMultipliers(*y_c, *y_d);
    }
    else {
      SmartPtr<Vector> y_c = IpData().curr_y_c()->MakeNew();
      SmartPtr<Vector> y_d = IpData().curr_y_d()->MakeNew();
      y_c->Set(0.0);
      y_d->Set(0.0);
      IpData().SetTrialEqMultipliers(*y_c, *y_d);
    }

    DBG_PRINT_VECTOR(2, "y_c", *IpData().curr_y_c());
    DBG_PRINT_VECTOR(2, "y_d", *IpData().curr_y_d());

    DBG_PRINT_VECTOR(2, "z_L", *IpData().curr_z_L());
    DBG_PRINT_VECTOR(2, "z_U", *IpData().curr_z_U());
    DBG_PRINT_VECTOR(2, "v_L", *IpData().curr_v_L());
    DBG_PRINT_VECTOR(2, "v_U", *IpData().curr_v_U());

    // upgrade the trial to the current point
    IpData().AcceptTrialPoint();

    return true;
  }

} // namespace Ipopt
