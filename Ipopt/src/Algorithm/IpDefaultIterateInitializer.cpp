// Copyright (C) 2004, 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter              IBM    2004-09-23

#include "IpDefaultIterateInitializer.hpp"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  DefaultIterateInitializer::DefaultIterateInitializer
  (const SmartPtr<EqMultiplierCalculator>& eq_mult_calculator,
   const SmartPtr<IterateInitializer>& warm_start_initializer,
   const SmartPtr<AugSystemSolver> aug_system_solver /*= NULL*/)
      :
      IterateInitializer(),
      eq_mult_calculator_(eq_mult_calculator),
      warm_start_initializer_(warm_start_initializer),
      aug_system_solver_(aug_system_solver)
  {}

  void DefaultIterateInitializer::RegisterOptions(SmartPtr<RegisteredOptions> reg_options)
  {
    reg_options->AddLowerBoundedNumberOption(
      "bound_push",
      "Desired minimum absolute distance from the initial point to bound.",
      0.0, true, 0.01,
      "Determines how much the initial point might have to "
      "be modified in order to be sufficiently inside "
      "the bounds (together with \"bound_frac\").  (This is kappa_1 in "
      "Section 3.6 of implementation paper.)");
    reg_options->AddBoundedNumberOption(
      "bound_frac",
      "Desired minimum relative distance from the initial point to bound.",
      0, true, 0.5, false, 0.01,
      "Determines how much the initial point might have to "
      "be modified in order to be sufficiently inside "
      "the bounds (together with \"bound_push\").  (This is kappa_2 in "
      "Section 3.6 of implementation paper.)");
    reg_options->AddLowerBoundedNumberOption(
      "slack_bound_push",
      "Desired minimum absolute distance from the initial slack to bound.",
      0.0, true, 0.01,
      "Determines how much the initial slack variables might have to "
      "be modified in order to be sufficiently inside "
      "the inequality bounds (together with \"slack_bound_frac\").  (This is kappa_1 in "
      "Section 3.6 of implementation paper.)");
    reg_options->AddBoundedNumberOption(
      "slack_bound_frac",
      "Desired minimum relative distance from the initial slack to bound.",
      0, true, 0.5, false, 0.01,
      "Determines how much the initial slack variables might have to "
      "be modified in order to be sufficiently inside "
      "the inequality bounds (together with \"slack_bound_push\").  (This is kappa_2 in "
      "Section 3.6 of implementation paper.)");
    reg_options->AddLowerBoundedNumberOption(
      "constr_mult_init_max",
      "Maximum allowed least-square guess of constraint multipliers.",
      0, false, 1e3,
      "Determines how large the initial least-square guesses of the constraint "
      "multipliers are allowed to be (in max-norm). If the guess is larger "
      "than this value, it is discarded and all constraint multipliers are "
      "set to zero.  This options is also used when initializing the "
      "restoration phase. By default, \"resto.constr_mult_init_max\" (the one "
      "used in RestoIterateInitializer) is set to zero.");
    reg_options->AddLowerBoundedNumberOption(
      "bound_mult_init_val",
      "Initial value for the bound multipliers.",
      0, true, 1.0,
      "All dual variables corresponding to bound constraints are "
      "initialized to this value.");
    reg_options->AddStringOption2(
      "bound_mult_init_method",
      "Initialization method for bound multipliers",
      "constant",
      "constant", "set all bound multipliers to the value of bound_mult_init_val",
      "mu-based", "initialize to mu_init/x_slack",
      "This option defines how the iterates for the bound multipliers are "
      "initialized.  If \"constant\" is chosen, then all bound multipliers "
      "are initialized to the value of \"bound_mult_init_val\".  If "
      "\"mu-based\" is chosen, the each value is initialized to the the value "
      "of \"mu_init\" divided by the corresponding slack variable.  This "
      "latter option might be useful if the starting point is close to the "
      "optimal solution.");
    reg_options->AddStringOption2(
      "least_square_init_primal",
      "Least square initialization of the primal variables", "no",
      "no", "take user-provided point",
      "yes", "overwrite user-provided point with least-square estimates",
      "If set to yes, Ipopt ignores the user provided point and solves a "
      "least square problem for the primal variables (x and s), to fit the "
      "linearized equality and inequality constraints.  This might be useful "
      "if the user doesn't know anything about the starting point, or for "
      "solving an LP or QP.");
    reg_options->AddStringOption2(
      "least_square_init_duals",
      "Least square initialization of all dual variables", "no",
      "no", "use bound_mult_init_val and least-square equality constraint multipliers",
      "yes", "overwrite user-provided point with least-square estimates",
      "If set to yes, Ipopt tries to compute least-square multipliers "
      "(considering ALL dual variables).  If successful, the bound "
      "multipliers are possibly corrected to be at least bound_mult_init_val. "
      "This might be useful "
      "if the user doesn't know anything about the starting point, or for "
      "solving an LP or QP.  This overwrites option "
      "\"bound_mult_init_method\".");
    reg_options->SetRegisteringCategory("Warm Start");
    reg_options->AddStringOption2(
      "warm_start_init_point",
      "Warm-start for initial point", "no",
      "no", "do not use the warm start initialization",
      "yes", "use the warm start initialization",
      "Indicates whether this optimization should use a warm start "
      "initialization, where values of primal and dual variables are "
      "given (e.g., from a previous optimization of a related problem.)");
  }

  bool DefaultIterateInitializer::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    // Check for the algorithm options
    options.GetNumericValue("bound_push", bound_push_, prefix);
    options.GetNumericValue("bound_frac", bound_frac_, prefix);
    if (!options.GetNumericValue("slack_bound_push", slack_bound_push_, prefix)) {
      slack_bound_push_ = bound_push_;
    }
    if (!options.GetNumericValue("slack_bound_frac", slack_bound_frac_, prefix)) {
      slack_bound_frac_ = bound_frac_;
    }
    options.GetNumericValue("constr_mult_init_max",
                            constr_mult_init_max_, prefix);
    options.GetNumericValue("bound_mult_init_val",
                            bound_mult_init_val_, prefix);
    options.GetBoolValue("warm_start_init_point",
                         warm_start_init_point_, prefix);
    options.GetBoolValue("least_square_init_primal",
                         least_square_init_primal_, prefix);
    ASSERT_EXCEPTION(!least_square_init_primal_ || IsValid(aug_system_solver_),
                     OPTION_INVALID,
                     "The least_square_init_primal can only be chosen if the DefaultInitializer object has an AugSystemSolver.\n");
    options.GetBoolValue("least_square_init_duals",
                         least_square_init_duals_, prefix);
    ASSERT_EXCEPTION(!least_square_init_duals_ || IsValid(aug_system_solver_),
                     OPTION_INVALID,
                     "The least_square_init_duals can only be chosen if the DefaultInitializer object has an AugSystemSolver.\n");
    int enum_int;
    options.GetEnumValue("bound_mult_init_method", enum_int, prefix);
    bound_mult_init_method_ = BoundMultInitMethod(enum_int);
    if (bound_mult_init_method_ == B_MU_BASED) {
      options.GetNumericValue("mu_init",  mu_init_, prefix);
    }

    bool retvalue = true;
    if (IsValid(eq_mult_calculator_)) {
      retvalue = eq_mult_calculator_->Initialize(Jnlst(), IpNLP(), IpData(),
                 IpCq(), options, prefix);
      if (!retvalue) {
        return retvalue;
      }
    }
    if (IsValid(warm_start_initializer_)) {
      retvalue =
        warm_start_initializer_->Initialize(Jnlst(), IpNLP(), IpData(),
                                            IpCq(), options, prefix);
    }
    return retvalue;
  }

  bool DefaultIterateInitializer::SetInitialIterates()
  {
    DBG_START_METH("DefaultIterateInitializer::SetInitialIterates",
                   dbg_verbosity);

    if (warm_start_init_point_) {
      DBG_ASSERT(IsValid(warm_start_initializer_));
      return warm_start_initializer_->SetInitialIterates();
    }

    // Get the starting values provided by the NLP and store them
    // in the ip_data current fields.  The following line only requests
    // intial values for the primal variables x, but later we might
    // make this more flexible based on user options.

    /////////////////////////////////////////////////////////////////////
    //                   Initialize primal variables                   //
    /////////////////////////////////////////////////////////////////////

    if (!IpData().InitializeDataStructures(IpNLP(), true, false, false,
                                           false, false)) {
      return false;
    }

    // get a container of the current point. We will modify parts of
    // this IteratesVector to set the trial point.
    SmartPtr<IteratesVector> iterates = IpData().curr()->MakeNewContainer();
    if (least_square_init_primal_) {
      // If least_square_init_primal, then we compute the least square x
      // and s that satisfy the linearized constraints, and push it into
      // bounds later.
      SmartPtr<Vector> x_ls = iterates->x()->MakeNew();
      SmartPtr<Vector> s_ls = iterates->s()->MakeNew();
      bool retval =CalculateLeastSquarePrimals(*x_ls, *s_ls);
      if (retval) {
        Jnlst().Printf(J_DETAILED, J_INITIALIZATION,
                       "Least square intial values for x and s computed.\n");
        x_ls->Print(Jnlst(), J_VECTOR, J_INITIALIZATION, "x_ls");
        s_ls->Print(Jnlst(), J_VECTOR, J_INITIALIZATION, "s_ls");
        iterates->Set_x(*x_ls);
        iterates->Set_s(*s_ls);
      }
      else {
        Jnlst().Printf(J_WARNING, J_INITIALIZATION,
                       "Least square initialization of x and s failed!\n");
      }
    }
    DBG_PRINT_VECTOR(2, "curr_x", *iterates->x());

    // Now we compute the initial values that the algorithm is going to
    // actually use.  We first store them in the trial fields in ip_data.

    // Push the x iterates sufficiently inside the bounds
    // Calculate any required shift in x0 and s0

    SmartPtr<const Vector> new_x;
    push_variables(Jnlst(), bound_push_, bound_frac_,
                   "x", *iterates->x(), new_x, *IpNLP().x_L(),
                   *IpNLP().x_U(), *IpNLP().Px_L(), *IpNLP().Px_U());

    iterates->Set_x(*new_x);
    IpData().set_trial(iterates);

    // Calculate the shift in s...
    SmartPtr<const Vector> s = IpCq().trial_d();
    DBG_PRINT_VECTOR(2, "s", *s);

    SmartPtr<const Vector> new_s;
    push_variables(Jnlst(), slack_bound_push_, slack_bound_frac_,
                   "s", *s, new_s, *IpNLP().d_L(),
                   *IpNLP().d_U(), *IpNLP().Pd_L(), *IpNLP().Pd_U());

    iterates = IpData().trial()->MakeNewContainer();
    iterates->Set_s(*new_s);

    /////////////////////////////////////////////////////////////////////
    //                   Initialize bound multipliers                  //
    /////////////////////////////////////////////////////////////////////

    // Initialize the bound multipliers to bound_mult_init_val.
    switch (bound_mult_init_method_) {
    case B_CONSTANT:
      iterates->create_new_z_L();
      iterates->create_new_z_U();
      iterates->create_new_v_L();
      iterates->create_new_v_U();
      iterates->z_L_NonConst()->Set(bound_mult_init_val_);
      iterates->z_U_NonConst()->Set(bound_mult_init_val_);
      iterates->v_L_NonConst()->Set(bound_mult_init_val_);
      iterates->v_U_NonConst()->Set(bound_mult_init_val_);
      IpData().set_trial(iterates);
      break;
    case B_MU_BASED:
      IpData().set_trial(iterates);
      iterates = IpData().trial()->MakeNewContainer();
      iterates->create_new_z_L();
      iterates->create_new_z_U();
      iterates->create_new_v_L();
      iterates->create_new_v_U();
      iterates->z_L_NonConst()->Set(mu_init_);
      iterates->z_U_NonConst()->Set(mu_init_);
      iterates->v_L_NonConst()->Set(mu_init_);
      iterates->v_U_NonConst()->Set(mu_init_);
      iterates->z_L_NonConst()->ElementWiseDivide(*IpCq().trial_slack_x_L());
      iterates->z_U_NonConst()->ElementWiseDivide(*IpCq().trial_slack_x_U());
      iterates->v_L_NonConst()->ElementWiseDivide(*IpCq().trial_slack_s_L());
      iterates->v_U_NonConst()->ElementWiseDivide(*IpCq().trial_slack_s_U());
      IpData().set_trial(iterates);
      break;
    default:
      THROW_EXCEPTION(OPTION_INVALID, "Invalid value of option bound_mult_init_method");
      break;
    }

    bool call_least_square_mults = true;
    if (least_square_init_duals_) {
      // We try to compute a least square estimate of all multiplers,
      // and if successful,we make sure they are sufficiently positive
      SmartPtr<Vector> zL_new = IpData().trial()->z_L()->MakeNew();
      SmartPtr<Vector> zU_new = IpData().trial()->z_U()->MakeNew();
      SmartPtr<Vector> vL_new = IpData().trial()->v_L()->MakeNew();
      SmartPtr<Vector> vU_new = IpData().trial()->v_U()->MakeNew();
      SmartPtr<Vector> yc_new = IpData().trial()->y_c()->MakeNew();
      SmartPtr<Vector> yd_new = IpData().trial()->y_d()->MakeNew();
      bool retval = CalculateLeastSquareDuals(*zL_new, *zU_new, *vL_new,
                                              *vU_new, *yc_new, *yd_new);
      if (retval) {
        // z_L etc are still at bound_mult_init_val_
        zL_new->ElementWiseMax(*IpData().trial()->z_L());
        zU_new->ElementWiseMax(*IpData().trial()->z_U());
        vL_new->ElementWiseMax(*IpData().trial()->v_L());
        vU_new->ElementWiseMax(*IpData().trial()->v_U());
        iterates = IpData().trial()->MakeNewContainer();
        iterates->Set_z_L(*zL_new);
        iterates->Set_z_U(*zU_new);
        iterates->Set_v_L(*vL_new);
        iterates->Set_v_U(*vU_new);
        iterates->Set_y_c(*yc_new);
        iterates->Set_y_d(*yd_new);
        IpData().set_trial(iterates);

        Jnlst().Printf(J_DETAILED, J_INITIALIZATION,
                       "Least square intial values for z_L, z_U,v_L, v_U, y_c, y_d computed.\n");
        zL_new->Print(Jnlst(), J_VECTOR, J_INITIALIZATION, "zL_new");
        zU_new->Print(Jnlst(), J_VECTOR, J_INITIALIZATION, "zU_new");
        vL_new->Print(Jnlst(), J_VECTOR, J_INITIALIZATION, "vL_new");
        vU_new->Print(Jnlst(), J_VECTOR, J_INITIALIZATION, "vU_new");
        yc_new->Print(Jnlst(), J_VECTOR, J_INITIALIZATION, "yc_new");
        yd_new->Print(Jnlst(), J_VECTOR, J_INITIALIZATION, "yd_new");
        call_least_square_mults = false;
      }
      else {
        Jnlst().Printf(J_WARNING, J_INITIALIZATION,
                       "Least square initialization of z_L, z_U,v_L, v_U, y_c, y_d failed!\n");
      }
    }
    if (call_least_square_mults) {
      /////////////////////////////////////////////////////////////////////
      //           Initialize equality constraint multipliers            //
      /////////////////////////////////////////////////////////////////////

      least_square_mults(Jnlst(), IpNLP(), IpData(), IpCq(),
                         eq_mult_calculator_, constr_mult_init_max_);
    }

    // upgrade the trial to the current point
    IpData().AcceptTrialPoint();

    return true;
  }

  bool
  DefaultIterateInitializer::CalculateLeastSquarePrimals(Vector& x_ls, Vector& s_ls)
  {
    DBG_START_METH("DefaultIterateInitializer::CalculateLeastSquarePrimals",
                   dbg_verbosity);
    SmartPtr<const SymMatrix> zeroW = IpNLP().uninitialized_h();
    DBG_PRINT_MATRIX(2, "zeroW", *zeroW);
    SmartPtr<const Matrix> J_c = IpCq().curr_jac_c();
    SmartPtr<const Matrix> J_d = IpCq().curr_jac_d();

    // Compute the right hand side
    SmartPtr<Vector> rhs_x = x_ls.MakeNew();
    rhs_x->Set(0.);
    SmartPtr<Vector> rhs_s = s_ls.MakeNew();
    rhs_s->Set(0.);

    SmartPtr<const Vector> rhs_c = IpCq().curr_c();
    SmartPtr<const Vector> rhs_d = IpCq().curr_d();

    SmartPtr<Vector> sol_c = rhs_c->MakeNew();
    SmartPtr<Vector> sol_d = rhs_d->MakeNew();

    DBG_PRINT_VECTOR(2, "rhs_x", *rhs_x);
    DBG_PRINT_VECTOR(2, "rhs_s", *rhs_s);
    DBG_PRINT_VECTOR(2, "rhs_c", *rhs_c);
    DBG_PRINT_VECTOR(2, "rhs_d", *rhs_d);

    enum ESymSolverStatus retval;
    Index numberOfEVals=rhs_c->Dim()+rhs_d->Dim();
    retval = aug_system_solver_->Solve(GetRawPtr(zeroW), 0.0, NULL, 1.0, NULL,
                                       1.0, GetRawPtr(J_c), NULL, 0.,
                                       GetRawPtr(J_d), NULL, 0., *rhs_x, *rhs_s,
                                       *rhs_c, *rhs_d, x_ls, s_ls, *sol_c, *sol_d,
                                       true, numberOfEVals);
    if (retval!=SYMSOLVER_SUCCESS) {
      return false;
    }
    x_ls.Scal(-1.);
    s_ls.Scal(-1.);

    DBG_PRINT_VECTOR(2, "sol_x", x_ls);
    DBG_PRINT_VECTOR(2, "sol_s", s_ls);
    DBG_PRINT_VECTOR(2, "sol_c", *sol_c);
    DBG_PRINT_VECTOR(2, "sol_d", *sol_d);

    return true;
  }

  bool
  DefaultIterateInitializer::
  CalculateLeastSquareDuals(Vector& zL_new, Vector& zU_new,
                            Vector& vL_new, Vector& vU_new,
                            Vector& yc_new, Vector& yd_new)
  {
    DBG_START_METH("DefaultIterateInitializer::CalculateLeastSquarePrimals",
                   dbg_verbosity);

    SmartPtr<const SymMatrix> zeroW = IpNLP().uninitialized_h();
    DBG_PRINT_MATRIX(2, "zeroW", *zeroW);
    SmartPtr<const Matrix> J_c = IpCq().curr_jac_c();
    SmartPtr<const Matrix> J_d = IpCq().curr_jac_d();

    // Compute the entries in Hessian diagonals
    SmartPtr<Vector> Dx = IpData().trial()->x()->MakeNew();
    SmartPtr<Vector> tmp = IpNLP().x_L()->MakeNew();
    tmp->Set(-1.);
    IpNLP().Px_L()->MultVector(1., *tmp, 0., *Dx);
    tmp = IpNLP().x_U()->MakeNew();
    tmp->Set(-1.);
    IpNLP().Px_U()->MultVector(1., *tmp, 1., *Dx);
    SmartPtr<Vector> Ds = IpData().trial()->s()->MakeNew();
    tmp = IpNLP().d_L()->MakeNew();
    tmp->Set(-1.);
    IpNLP().Pd_L()->MultVector(1., *tmp, 0., *Ds);
    tmp = IpNLP().d_U()->MakeNew();
    tmp->Set(-1.);
    IpNLP().Pd_U()->MultVector(1., *tmp, 1., *Ds);

    // Get the right hand side
    SmartPtr<const Vector> rhs_x = IpCq().trial_grad_f();
    SmartPtr<Vector> rhs_s = Ds->MakeNew();
    rhs_s->Set(0.);
    SmartPtr<Vector> rhs_c = yc_new.MakeNew();
    rhs_c->Set(0.);
    SmartPtr<Vector> rhs_d = yd_new.MakeNew();
    rhs_d->Set(0.);

    // Space for the solution
    SmartPtr<Vector> sol_x = rhs_x->MakeNew();
    SmartPtr<Vector> sol_s = rhs_s->MakeNew();

    DBG_PRINT_VECTOR(2, "rhs_x", *rhs_x);
    DBG_PRINT_VECTOR(2, "rhs_s", *rhs_s);
    DBG_PRINT_VECTOR(2, "rhs_c", *rhs_c);
    DBG_PRINT_VECTOR(2, "rhs_d", *rhs_d);

    enum ESymSolverStatus retval;
    Index numberOfEVals=rhs_x->Dim()+rhs_s->Dim();
    retval =
      aug_system_solver_->Solve(GetRawPtr(zeroW), 0.0, GetRawPtr(Dx), 0.0, GetRawPtr(Ds),
                                0.0, GetRawPtr(J_c), NULL, 0.,
                                GetRawPtr(J_d), NULL, 0., *rhs_x, *rhs_s,
                                *rhs_c, *rhs_d, *sol_x, *sol_s,
                                yc_new, yd_new, true, numberOfEVals);
    if (retval!=SYMSOLVER_SUCCESS) {
      return false;
    }
    DBG_PRINT_VECTOR(2, "sol_x", *sol_x);
    DBG_PRINT_VECTOR(2, "sol_s", *sol_s);
    DBG_PRINT_VECTOR(2, "sol_c", yc_new);
    DBG_PRINT_VECTOR(2, "sol_d", yd_new);

    // Get the output right
    yc_new.Scal(-1.0);
    yd_new.Scal(-1.0);
    IpNLP().Px_L()->TransMultVector(-1., *sol_x, 0., zL_new);
    IpNLP().Px_U()->TransMultVector(1., *sol_x, 0., zU_new);
    IpNLP().Pd_L()->TransMultVector(-1., *sol_s, 0., vL_new);
    IpNLP().Pd_U()->TransMultVector(1., *sol_s, 0., vU_new);

    return true;
  }

  void DefaultIterateInitializer::push_variables(
    const Journalist& jnlst,
    Number bound_push,
    Number bound_frac,
    std::string name,
    const Vector& orig_x,
    SmartPtr<const Vector>& new_x,
    const Vector& x_L,
    const Vector& x_U,
    const Matrix& Px_L,
    const Matrix& Px_U)
  {
    DBG_START_FUN("DefaultIterateInitializer::push_variables",
                  dbg_verbosity);

    SmartPtr<const Vector> my_orig_x = &orig_x;
    // ToDo: Make this more efficient...?

    // To avoid round-off error, move variables first at the bounds
    if (bound_push>0. || bound_frac>0.) {
      push_variables(jnlst, 0., 0., name, orig_x, new_x, x_L, x_U, Px_L, Px_U);
      my_orig_x = new_x;
    }

    DBG_PRINT_VECTOR(2,"orig_x", *my_orig_x);
    DBG_PRINT_MATRIX(2,"Px_L", Px_L);
    DBG_PRINT_VECTOR(2, "x_L", x_L);
    DBG_PRINT_MATRIX(2,"Px_U", Px_U);
    DBG_PRINT_VECTOR(2, "x_U", x_U);

    SmartPtr<Vector> tmp_l = x_L.MakeNew();
    SmartPtr<Vector> tmp_u = x_U.MakeNew();

    const double dbl_min = std::numeric_limits<double>::min();
    const double tiny_double = 100.0*dbl_min;

    // Calculate any required shift in x0 and s0
    SmartPtr<Vector> tmp = my_orig_x->MakeNew();
    SmartPtr<Vector> tiny_l = x_L.MakeNew();
    tiny_l->Set(tiny_double);

    SmartPtr<Vector> q_l = x_L.MakeNew();
    SmartPtr<Vector> p_l = x_L.MakeNew();
    SmartPtr<Vector> delta_x = my_orig_x->MakeNew();

    SmartPtr<Vector> zero_l = x_L.MakeNew();
    zero_l->Set(0.0);
    SmartPtr<Vector> zero_u = x_U.MakeNew();
    zero_u->Set(0.0);

    if (bound_frac>0.) {
      DBG_ASSERT(bound_push>0.);

      // Calculate p_l
      Px_L.MultVector(1.0, x_L, 0.0, *tmp);
      Px_U.TransMultVector(1.0, *tmp, 0.0, *tmp_u);
      tmp_u->AddOneVector(1., x_U, -1.);
      Px_U.MultVector(1.0, *tmp_u, 0.0, *tmp);
      Px_L.TransMultVector(1.0, *tmp, 0.0, *q_l);
      q_l->AddOneVector(-1.0, *tiny_l, bound_frac);
      DBG_PRINT_VECTOR(2, "q_l", *q_l);
      // At this point, q_l is
      // bound_frac * Px_L^T Px_U(x_U - Px_U^T Px_L x_L)  -  tiny_double
      // i.e., it is bound_frac*(x_U - x_L) for those components that have
      //          two bounds
      //       and - tiny_double for those that have only one or no bounds

      tmp_l->Set(bound_push);
      p_l->AddOneVector(bound_push, x_L, 0.);
      p_l->ElementWiseAbs();
      p_l->ElementWiseMax(*tmp_l);
      // now p_l is bound_push * max(|x_L|,1)

      q_l->ElementWiseReciprocal();
      p_l->ElementWiseReciprocal();

      p_l->ElementWiseMax(*q_l);
      p_l->ElementWiseReciprocal();
      //    p_l->Axpy(1.0, *tiny_l);  we shouldn't need this

      // At this point, p_l is
      //  min(bound_push * max(|x_L|,1), bound_frac*(x_U-x_L)) for components
      //                                                       with two bounds
      //  bound_push * max(|x_L|,1)                            otherwise
      // This is the margin we want to the lower bound
      DBG_PRINT_VECTOR(1, "p_l", *p_l);

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
      q_u->AddOneVector(-1.0, *tiny_u, bound_frac);
      DBG_PRINT_VECTOR(2,"q_u",*q_u);
      // q_u is now the same as q_l above, but of the same dimension as x_L

      tmp_u->Set(bound_push);
      p_u->Copy(x_U);
      p_u->AddOneVector(bound_push, x_U, 0.);
      p_u->ElementWiseAbs();
      p_u->ElementWiseMax(*tmp_u);
      DBG_PRINT_VECTOR(2,"p_u",*p_u);

      q_u->ElementWiseReciprocal();
      p_u->ElementWiseReciprocal();

      p_u->ElementWiseMax(*q_u);
      p_u->ElementWiseReciprocal();
      p_u->Axpy(1.0, *tiny_u);
      // At this point, p_l is
      //  min(bound_push * max(|x_U|,1), bound_frac*(x_U-x_L)) for components
      //                                                       with two bounds
      //  bound_push * max(|x_U|,1)                            otherwise
      // This is the margin we want to the upper bound
      DBG_PRINT_VECTOR(2,"actual_p_u",*p_u);

      // Calculate the new x

      Px_L.TransMultVector(-1.0, *my_orig_x, 0.0, *tmp_l);
      DBG_PRINT_VECTOR(1, "tmp_l1", *tmp_l);
      tmp_l->AddTwoVectors(1.0, x_L, 1.0, *p_l, 1.);
      tmp_l->ElementWiseMax(*zero_l);
      DBG_PRINT_VECTOR(1, "tmp_l2", *tmp_l);
      // tmp_l is now max(x_L + p_l - x, 0), i.e., the amount by how
      // much need to correct the variable

      Px_U.TransMultVector(1.0, *my_orig_x, 0.0, *tmp_u);
      tmp_u->AddTwoVectors(-1.0, x_U, 1.0, *p_u, 1.);
      tmp_u->ElementWiseMax(*zero_u);
      // tmp_u is now max(x - (x_U-p_u), 0), i.e., the amount by how
      // much need to correct the variable
    }
    else {
      DBG_ASSERT(bound_push == 0.);
      tmp_l = x_L.MakeNewCopy();
      DBG_PRINT_VECTOR(1, "tmp_l33", *tmp_l);
      Px_L.TransMultVector(-1.0, *my_orig_x, 1.0, *tmp_l);
      tmp_l->ElementWiseMax(*zero_l);
      DBG_PRINT_VECTOR(1, "tmp_l3", *tmp_l);

      tmp_u = x_U.MakeNewCopy();
      Px_U.TransMultVector(1.0, *my_orig_x, -1.0, *tmp_u);
      tmp_u->ElementWiseMax(*zero_u);
    }

    Number nrm_l = tmp_l->Amax();
    if (nrm_l>0.) {
      Px_L.MultVector(1.0, *tmp_l, 0.0, *delta_x);
    }
    else {
      delta_x->Set(0.);
    }
    Number nrm_u = tmp_u->Amax();
    if (nrm_u>0.) {
      Px_U.MultVector(-1.0, *tmp_u, 1.0, *delta_x);
    }

    if (nrm_l > 0 || nrm_u > 0) {
      delta_x->Axpy(1.0, *my_orig_x);
      new_x = ConstPtr(delta_x);
      if (bound_push > 0.) {
        jnlst.Printf(J_DETAILED, J_INITIALIZATION,
                     "Moved initial values of %s sufficiently inside the bounds.\n", name.c_str());
        my_orig_x->Print(jnlst, J_VECTOR, J_INITIALIZATION, "original vars");
        new_x->Print(jnlst, J_VECTOR, J_INITIALIZATION, "new vars");
      }
    }
    else {
      new_x = my_orig_x;
      if (bound_push > 0.) {
        jnlst.Printf(J_DETAILED, J_INITIALIZATION,
                     "Initial values of %s sufficiently inside the bounds.\n", name.c_str());
      }
    }
  }

  void DefaultIterateInitializer::least_square_mults(
    const Journalist& jnlst,
    IpoptNLP& ip_nlp,
    IpoptData& ip_data,
    IpoptCalculatedQuantities& ip_cq,
    const SmartPtr<EqMultiplierCalculator>& eq_mult_calculator,
    Number constr_mult_init_max)
  {
    DBG_START_FUN("DefaultIterateInitializer::least_square_mults",
                  dbg_verbosity);

    SmartPtr<IteratesVector> iterates = ip_data.trial()->MakeNewContainer();
    iterates->create_new_y_c();
    iterates->create_new_y_d();

    if (iterates->y_c_NonConst()->Dim()==iterates->x()->Dim()) {
      // This problem is square, we just set the multipliers to zero
      iterates->y_c_NonConst()->Set(0.0);
      iterates->y_d_NonConst()->Set(0.0);
      ip_data.Append_info_string("s");
    }
    else if (IsValid(eq_mult_calculator) && constr_mult_init_max>0. &&
             iterates->y_c_NonConst()->Dim()+iterates->y_d_NonConst()->Dim()>0) {
      // First move all the trial data into the current fields, since
      // those values are needed to compute the initial values for
      // the multipliers
      ip_data.CopyTrialToCurrent();
      SmartPtr<Vector> y_c = iterates->y_c_NonConst();
      SmartPtr<Vector> y_d = iterates->y_d_NonConst();
      bool retval = eq_mult_calculator->CalculateMultipliers(*y_c, *y_d);
      if (!retval) {
        y_c->Set(0.0);
        y_d->Set(0.0);
      }
      else {
        /*
        {
          ip_data.set_trial(iterates);
          printf("grad_x = %e grad_s = %e y_c = %e y_d = %e\n",
          ip_cq.trial_grad_lag_x()->Amax(),
          ip_cq.trial_grad_lag_s()->Amax(),
          y_c->Amax(),
          y_d->Amax());
          iterates = ip_data.trial()->MakeNewContainer();
        }
        */
        jnlst.Printf(J_DETAILED, J_INITIALIZATION,
                     "Least square estimates max(y_c) = %e, max(y_d) = %e\n",
                     y_c->Amax(), y_d->Amax());
        Number yinitnrm = Max(y_c->Amax(), y_d->Amax());
        if (yinitnrm > constr_mult_init_max) {
          y_c->Set(0.0);
          y_d->Set(0.0);
        }
        else {
          ip_data.Append_info_string("y");
        }
      }
    }
    else {
      iterates->y_c_NonConst()->Set(0.0);
      iterates->y_d_NonConst()->Set(0.0);
    }
    ip_data.set_trial(iterates);

    DBG_PRINT_VECTOR(2, "y_c", *ip_data.trial()->y_c());
    DBG_PRINT_VECTOR(2, "y_d", *ip_data.trial()->y_d());
  }

} // namespace Ipopt
