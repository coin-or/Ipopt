// Copyright (C) 2004, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter              IBM    2004-10-12

#include "IpRestoIterateInitializer.hpp"
#include "IpRestoIpoptNLP.hpp"
#include "IpDefaultIterateInitializer.hpp"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  RestoIterateInitializer::RestoIterateInitializer
  (const SmartPtr<EqMultiplierCalculator>& resto_eq_mult_calculator)
      :
      IterateInitializer(),
      resto_eq_mult_calculator_(resto_eq_mult_calculator)
  {}

  bool RestoIterateInitializer::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    if (!options.GetNumericValue("constr_mult_init_max", constr_mult_init_max_, prefix)) {
      // By default, we want to set this to zero. Seems to work better
      // as initialization for the restoration phase
      constr_mult_init_max_ = 0.;
    }

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
      static_cast<RestoIpoptNLP*> (&IpNLP());
    SmartPtr<IpoptNLP> orig_ip_nlp =
      static_cast<IpoptNLP*> (&resto_ip_nlp->OrigIpNLP());
    SmartPtr<IpoptData> orig_ip_data =
      static_cast<IpoptData*> (&resto_ip_nlp->OrigIpData());
    SmartPtr<IpoptCalculatedQuantities> orig_ip_cq =
      static_cast<IpoptCalculatedQuantities*> (&resto_ip_nlp->OrigIpCq());

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
                                      false, false);

    SmartPtr<Vector> new_x = IpData().curr()->x()->MakeNew();
    SmartPtr<CompoundVector> Cnew_x =
      static_cast<CompoundVector*> (GetRawPtr(new_x));
    DBG_ASSERT(dynamic_cast<CompoundVector*> (GetRawPtr(new_x)));

    // Set the trial x variables from the original NLP
    Cnew_x->GetCompNonConst(0)->Copy(*orig_ip_data->curr()->x());

    // Compute the initial values for the n and p variables for the
    // equality constraints
    Number rho = resto_ip_nlp->Rho();
    DBG_PRINT((1,"rho = %e\n", rho));
    SmartPtr<Vector> nc = Cnew_x->GetCompNonConst(1);
    SmartPtr<Vector> pc = Cnew_x->GetCompNonConst(2);
    SmartPtr<const Vector> cvec = orig_ip_cq->curr_c();
    DBG_PRINT_VECTOR(2, "cvec", *cvec);
    SmartPtr<Vector> a = nc->MakeNew();
    SmartPtr<Vector> b = nc->MakeNew();
    a->Set(resto_mu/(2.*rho));
    a->Axpy(-0.5, *cvec);
    b->Copy(*cvec);
    b->Scal(resto_mu/(2.*rho));
    DBG_PRINT_VECTOR(2, "a", *a);
    DBG_PRINT_VECTOR(2, "b", *b);
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
    SmartPtr<Vector> new_s = IpData().curr()->s()->MakeNew();
    SmartPtr<CompoundVector> Cnew_s =
      static_cast<CompoundVector*> (GetRawPtr(new_s));
    DBG_ASSERT(dynamic_cast<CompoundVector*> (GetRawPtr(new_s)));
    DBG_ASSERT(Cnew_s->NComps() == 1);

    // Set the trial s variables from the original NLP
    Cnew_s->GetCompNonConst(0)->Copy(*orig_ip_data->curr()->s());

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
      static_cast<CompoundVector*> (GetRawPtr(new_z_L));
    DBG_ASSERT(dynamic_cast<CompoundVector*> (GetRawPtr(new_z_L)));
    SmartPtr<Vector> new_z_U = IpData().curr()->z_U()->MakeNew();
    SmartPtr<CompoundVector> Cnew_z_U =
      static_cast<CompoundVector*> (GetRawPtr(new_z_U));
    DBG_ASSERT(dynamic_cast<CompoundVector*> (GetRawPtr(new_z_U)));
    DBG_ASSERT(Cnew_z_U->NComps() == 1);
    SmartPtr<Vector> new_v_L = IpData().curr()->v_L()->MakeNew();
    SmartPtr<CompoundVector> Cnew_v_L =
      static_cast<CompoundVector*> (GetRawPtr(new_v_L));
    DBG_ASSERT(dynamic_cast<CompoundVector*> (GetRawPtr(new_v_L)));
    DBG_ASSERT(Cnew_v_L->NComps() == 1);
    SmartPtr<Vector> new_v_U = IpData().curr()->v_U()->MakeNew();
    SmartPtr<CompoundVector> Cnew_v_U =
      static_cast<CompoundVector*> (GetRawPtr(new_v_U));
    DBG_ASSERT(dynamic_cast<CompoundVector*> (GetRawPtr(new_v_U)));
    DBG_ASSERT(Cnew_v_U->NComps() == 1);

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
    SmartPtr<Vector> Cnew_z_U0 = Cnew_z_U->GetCompNonConst(0);
    Cnew_z_U0->Set(rho);
    Cnew_z_U0->ElementWiseMin(*orig_z_U);
    SmartPtr<Vector> Cnew_v_L0 = Cnew_v_L->GetCompNonConst(0);
    Cnew_v_L0->Set(rho);
    Cnew_v_L0->ElementWiseMin(*orig_v_L);
    SmartPtr<Vector> Cnew_v_U0 = Cnew_v_U->GetCompNonConst(0);
    Cnew_v_U0->Set(rho);
    Cnew_v_U0->ElementWiseMin(*orig_v_U);

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

    DefaultIterateInitializer::least_square_mults(
      Jnlst(), IpNLP(), IpData(), IpCq(),
      resto_eq_mult_calculator_, constr_mult_init_max_);

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
