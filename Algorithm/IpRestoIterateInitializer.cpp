// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter              IBM    2004-10-12

#include "IpRestoIterateInitializer.hpp"
#include "IpRestoIpoptNLP.hpp"

namespace Ipopt
{
  DBG_SET_VERBOSITY(0);

  DefineIpoptType(RestoIterateInitializer);

  RestoIterateInitializer::RestoIterateInitializer
  (const SmartPtr<EqMultiplierCalculator>& resto_eq_mult_calculator)
      :
      IterateInitializer(),
      resto_eq_mult_calculator_(resto_eq_mult_calculator)
  {}

  void RestoIterateInitializer::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddLowerBoundedNumberOption("laminitmax", "maximum value for the initial lambda's",
                                          0.0, false, 1e3);
  }

  bool RestoIterateInitializer::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    options.GetNumericValue("laminitmax", laminitmax_, prefix);

    bool retvalue = true;
    if (IsValid(resto_eq_mult_calculator_)) {
      retvalue = resto_eq_mult_calculator_->Initialize(Jnlst(),
                 IpNLP(), IpData(),
                 IpCq(), options, prefix);
    }
    return retvalue;
  }

  bool RestoIterateInitializer::SetInitialIterates()
  {
    DBG_START_METH("RestoIterateInitializer::SetInitialIterates",
                   dbg_verbosity);

    // Get a grip on the restoration phase NLP and obtain the pointers
    // to the original NLP data
    SmartPtr<RestoIpoptNLP> resto_ip_nlp =
      dynamic_cast<RestoIpoptNLP*> (&IpNLP());
    DBG_ASSERT(IsValid(resto_ip_nlp));
    SmartPtr<IpoptNLP> orig_ip_nlp =
      dynamic_cast<IpoptNLP*> (&resto_ip_nlp->OrigIpNLP());
    DBG_ASSERT(IsValid(orig_ip_nlp));
    SmartPtr<IpoptData> orig_ip_data =
      dynamic_cast<IpoptData*> (&resto_ip_nlp->OrigIpData());
    DBG_ASSERT(IsValid(orig_ip_data));
    SmartPtr<IpoptCalculatedQuantities> orig_ip_cq =
      dynamic_cast<IpoptCalculatedQuantities*> (&resto_ip_nlp->OrigIpCq());
    DBG_ASSERT(IsValid(orig_ip_cq));

    // Set the value of the barrier parameter
    Number resto_mu;
    resto_mu = Max(orig_ip_data->curr_mu(),
                   orig_ip_cq->curr_c()->Amax(),
                   orig_ip_cq->curr_d_minus_s()->Amax());
    IpData().Set_mu(resto_mu);
    Jnlst().Printf(J_DETAILED, J_INITIALIZATION,
                   "Initial barrier parameter resto_mu = %e\n", resto_mu);

    /////////////////////////////////////////////////////////////////////
    //                   Initialize primal varialbes                   //
    /////////////////////////////////////////////////////////////////////

    // initialize the data structures in the restoration phase NLP
    IpData().InitializeDataStructures(IpNLP(), false, false, false,
                                      false, false, false, false);

    SmartPtr<Vector> new_x = IpData().curr()->x()->MakeNew();
    SmartPtr<CompoundVector> Cnew_x =
      dynamic_cast<CompoundVector*> (GetRawPtr(new_x));

    // Set the trial x variables from the original NLP
    Cnew_x->GetCompNonConst(0)->Copy(*orig_ip_data->curr()->x());

    // Compute the initial values for the n and p variables for the
    // equality constraints
    Number rho = resto_ip_nlp->Rho();
    SmartPtr<Vector> nc = Cnew_x->GetCompNonConst(1);
    SmartPtr<Vector> pc = Cnew_x->GetCompNonConst(2);
    SmartPtr<const Vector> cvec = orig_ip_cq->curr_c();
    SmartPtr<Vector> a = nc->MakeNew();
    SmartPtr<Vector> b = nc->MakeNew();
    a->Set(resto_mu/(2.*rho));
    a->Axpy(-0.5, *cvec);
    b->Copy(*cvec);
    b->Scal(resto_mu/(2.*rho));
    solve_quadratic(*a, *b, *nc);
    pc->Copy(*cvec);
    pc->Axpy(1., *nc);
    DBG_PRINT_VECTOR(2, "nc", *nc);
    DBG_PRINT_VECTOR(2, "pc", *pc);

    // initial values for the n and p variables for the inequality
    // constraints
    SmartPtr<Vector> nd = Cnew_x->GetCompNonConst(3);
    SmartPtr<Vector> pd = Cnew_x->GetCompNonConst(4);
    cvec = orig_ip_cq->curr_d_minus_s();
    a = nd->MakeNew();
    b = nd->MakeNew();
    a->Set(resto_mu/(2.*rho));
    a->Axpy(-0.5, *cvec);
    b->Copy(*cvec);
    b->Scal(resto_mu/(2.*rho));
    solve_quadratic(*a, *b, *nd);
    pd->Copy(*cvec);
    pd->Axpy(1., *nd);
    DBG_PRINT_VECTOR(2, "nd", *nd);
    DBG_PRINT_VECTOR(2, "pd", *pd);

    // Leave the slacks unchanged
    SmartPtr<const Vector> new_s = orig_ip_data->curr()->s();

#ifdef orig
    // The initial values for the inequality n and p variables are
    // trickier, since only some constraints have both n and p
    // variables.  If they don't have both, the slack takes the role
    // of the missing one.
    SmartPtr<const Vector> d_L = orig_ip_nlp->d_L();
    SmartPtr<const Matrix> Pd_L = orig_ip_nlp->Pd_L();
    SmartPtr<const Vector> d_U = orig_ip_nlp->d_U();
    SmartPtr<const Matrix> Pd_U = orig_ip_nlp->Pd_U();
    // compute indicator vectors, that are 1. for those entries in d
    // that have only lower, only upper, or both bounds, respectively
    SmartPtr<Vector> ind_only_L = orig_ip_data->curr()->y_d()->MakeNew();
    SmartPtr<Vector> ind_only_U = orig_ip_data->curr()->y_d()->MakeNew();
    SmartPtr<Vector> ind_both = orig_ip_data->curr()->y_d()->MakeNew();
    SmartPtr<Vector> tmp = d_U->MakeNew();
    tmp->Set(1.);
    Pd_U->MultVector(1., *tmp, 0., *ind_only_U);
    tmp = d_L->MakeNew();
    tmp->Set(1.);
    Pd_L->MultVector(1., *tmp, 0., *ind_only_L);
    Pd_L->TransMultVector(1., *ind_only_U, 0., *tmp);
    Pd_L->MultVector(1., *tmp, 0., *ind_both);
    ind_only_L->Axpy(-1., *ind_both);
    ind_only_U->Axpy(-1., *ind_both);
    tmp = NULL; // free memory
    DBG_PRINT_VECTOR(2, "ind_only_L", *ind_only_L);
    DBG_PRINT_VECTOR(2, "ind_only_U", *ind_only_U);
    DBG_PRINT_VECTOR(2, "ind_both", *ind_both);

    // We now compute cvec to be
    // a) d - s    for entries with both bounds
    // b) d - d_L  for entries with only a lower bound
    // c) d - d_U  for entries with only an upper bound
    SmartPtr<Vector> cvec = ind_both->MakeNew();
    cvec->Copy(*orig_ip_cq->curr_d());
    DBG_PRINT_VECTOR(2,"orig d",*cvec);

    // a)
    SmartPtr<Vector> tmpfull = ind_both->MakeNew();
    tmpfull->Copy(*orig_ip_data->curr()->s());
    DBG_PRINT_VECTOR(1,"curr_s",*tmpfull);
    tmpfull->ElementWiseMultiply(*ind_both);
    DBG_PRINT_VECTOR(2,"s both",*tmpfull);
    cvec->Axpy(-1., *tmpfull);
    // b)
    Pd_L->MultVector(1., *d_L, 0., *tmpfull);
    tmpfull->ElementWiseMultiply(*ind_only_L);
    DBG_PRINT_VECTOR(2,"d_L only lower",*tmpfull);
    cvec->Axpy(-1., *tmpfull);
    // c)
    Pd_U->MultVector(1., *d_U, 0., *tmpfull);
    tmpfull->ElementWiseMultiply(*ind_only_U);
    DBG_PRINT_VECTOR(2,"d_U only upper",*tmpfull);
    cvec->Axpy(-1., *tmpfull);

    // now set up coefficients for the quadratic equation and solve it
    b = cvec->MakeNew();
    b->Set(1.);
    b->Axpy(1., *ind_both);
    b->ElementWiseReciprocal();
    b->Scal(resto_mu/rho);
    DBG_PRINT_VECTOR(1, "b2", *b);
    // now, b is resto_mu/rho for all entries with only one bound, and
    // resto_mu/(2rho) for all entries with two bounds
    a = cvec->MakeNew();
    a->Copy(*cvec);
    a->Scal(-0.5);
    a->Axpy(1., *b);
    b->ElementWiseMultiply(*cvec);
    solve_quadratic(*a, *b, *tmpfull);
    SmartPtr<Vector> tmpfull2 = tmpfull->MakeNew();
    tmpfull2->Copy(*cvec);
    tmpfull2->Axpy(1., *tmpfull);
    SmartPtr<Vector> nd = Cnew_x->GetCompNonConst(3);
    SmartPtr<Vector> pd = Cnew_x->GetCompNonConst(4);
    Pd_L->TransMultVector(1., *tmpfull, 0., *nd);
    Pd_U->TransMultVector(1., *tmpfull2, 0., *pd);
    DBG_PRINT_VECTOR(2,"nfull",*tmpfull);
    DBG_PRINT_VECTOR(2,"pfull",*tmpfull2);
    DBG_PRINT_VECTOR(2,"new nd",*nd);
    DBG_PRINT_VECTOR(2,"new pd",*pd);

    // the new values for the slacks is takes as
    // a) keep for those entries that have two bounds
    // b) set to pfull + d_L for those with only lower bounds
    // c) set to d_u - nfull for those with only upper bounds
    SmartPtr<Vector> new_s = IpData().curr()->s()->MakeNew();
    // a)
    new_s->Copy(*orig_ip_data->curr()->s());
    DBG_PRINT_VECTOR(2, "curr_s", *new_s);
    new_s->ElementWiseMultiply(*ind_both);
    // b)
    Pd_L->MultVector(1., *d_L, 1., *tmpfull2);
    tmpfull2->ElementWiseMultiply(*ind_only_L);
    new_s->Axpy(1., *tmpfull2);
    // c)
    Pd_U->MultVector(1., *d_U, -1., *tmpfull);
    tmpfull->ElementWiseMultiply(*ind_only_U);
    new_s->Axpy(1., *tmpfull);
#endif

    // Now set the primal trial variables
    DBG_PRINT_VECTOR(2,"new_s",*new_s);
    DBG_PRINT_VECTOR(2,"new_x",*new_x);
    SmartPtr<IteratesVector> trial = IpData().curr()->MakeNewContainer();
    trial->Set_primal(*new_x, *new_s);
    IpData().set_trial(trial);

    DBG_PRINT_VECTOR(2, "resto_c", *IpCq().trial_c());
    DBG_PRINT_VECTOR(2, "resto_d_minus_s", *IpCq().trial_d_minus_s());

    /////////////////////////////////////////////////////////////////////
    //                   Initialize bound multipliers                  //
    /////////////////////////////////////////////////////////////////////

    SmartPtr<Vector> new_z_L = IpData().curr()->z_L()->MakeNew();
    SmartPtr<CompoundVector> Cnew_z_L =
      dynamic_cast<CompoundVector*> (GetRawPtr(new_z_L));
    DBG_ASSERT(IsValid(Cnew_z_L));
    SmartPtr<Vector> new_z_U = IpData().curr()->z_U()->MakeNew();
    SmartPtr<Vector> new_v_L = IpData().curr()->v_L()->MakeNew();
    SmartPtr<Vector> new_v_U = IpData().curr()->v_U()->MakeNew();

    // multipliers for the original bounds are
    SmartPtr<const Vector> orig_z_L = orig_ip_data->curr()->z_L();
    SmartPtr<const Vector> orig_z_U = orig_ip_data->curr()->z_U();
    SmartPtr<const Vector> orig_v_L = orig_ip_data->curr()->v_L();
    SmartPtr<const Vector> orig_v_U = orig_ip_data->curr()->v_U();

    // Set the new multipliers to the min of the penalty parameter Rho
    // and their current value
    SmartPtr<Vector> Cnew_z_L0 = Cnew_z_L->GetCompNonConst(0);
    Cnew_z_L0->Set(rho);
    Cnew_z_L0->ElementWiseMin(*orig_z_L);
    new_z_U->Set(rho);
    new_z_U->ElementWiseMin(*orig_z_U);
    new_v_L->Set(rho);
    new_v_L->ElementWiseMin(*orig_v_L);
    new_v_U->Set(rho);
    new_v_U->ElementWiseMin(*orig_v_U);

    // Set the multipliers for the p and n bounds to the "primal" multipliers
    SmartPtr<Vector> Cnew_z_L1 = Cnew_z_L->GetCompNonConst(1);
    Cnew_z_L1->Set(resto_mu);
    Cnew_z_L1->ElementWiseDivide(*nc);
    SmartPtr<Vector> Cnew_z_L2 = Cnew_z_L->GetCompNonConst(2);
    Cnew_z_L2->Set(resto_mu);
    Cnew_z_L2->ElementWiseDivide(*pc);
    SmartPtr<Vector> Cnew_z_L3 = Cnew_z_L->GetCompNonConst(3);
    Cnew_z_L3->Set(resto_mu);
    Cnew_z_L3->ElementWiseDivide(*nd);
    SmartPtr<Vector> Cnew_z_L4 = Cnew_z_L->GetCompNonConst(4);
    Cnew_z_L4->Set(resto_mu);
    Cnew_z_L4->ElementWiseDivide(*pd);

    // Set those initial values to be the trial values in Data
    trial = IpData().trial()->MakeNewContainer();
    trial->Set_bound_mult(*new_z_L, *new_z_U, *new_v_L, *new_v_U);
    IpData().set_trial(trial);

    /////////////////////////////////////////////////////////////////////
    //           Initialize equality constraint multipliers            //
    /////////////////////////////////////////////////////////////////////

    if (IsValid(resto_eq_mult_calculator_) && laminitmax_>0.) {
      // First move all the trial data into the current fields, since
      // those values are needed to compute the initial values for
      // the multipliers
      IpData().CopyTrialToCurrent();
      trial = IpData().trial()->MakeNewContainer();
      SmartPtr<Vector> y_c = IpData().curr()->y_c()->MakeNew();
      SmartPtr<Vector> y_d = IpData().curr()->y_d()->MakeNew();
      bool retval = resto_eq_mult_calculator_->CalculateMultipliers(*y_c, *y_d);
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
      trial->Set_eq_mult(*y_c, *y_d);
    }
    else {
      SmartPtr<Vector> y_c = IpData().curr()->y_c()->MakeNew();
      SmartPtr<Vector> y_d = IpData().curr()->y_d()->MakeNew();
      y_c->Set(0.0);
      y_d->Set(0.0);
      trial->Set_eq_mult(*y_c, *y_d);
    }
    // update the new multipliers
    IpData().set_trial(trial);

    // upgrade the trial to the current point
    IpData().AcceptTrialPoint();

    DBG_PRINT_VECTOR(2, "y_c", *IpData().curr()->y_c());
    DBG_PRINT_VECTOR(2, "y_d", *IpData().curr()->y_d());

    DBG_PRINT_VECTOR(2, "z_L", *IpData().curr()->z_L());
    DBG_PRINT_VECTOR(2, "z_U", *IpData().curr()->z_U());
    DBG_PRINT_VECTOR(2, "v_L", *IpData().curr()->v_L());
    DBG_PRINT_VECTOR(2, "v_U", *IpData().curr()->v_U());

    return true;
  }

  void
  RestoIterateInitializer::solve_quadratic(const Vector& a,
      const Vector& b,
      Vector& v)
  {
    v.Copy(a);
    v.ElementWiseMultiply(a);

    v.Axpy(1., b);
    v.ElementWiseSqrt();

    v.Axpy(1., a);
  }

} // namespace Ipopt
