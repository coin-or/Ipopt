// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpIpoptAlg.hpp"
#include "IpJournalist.hpp"

namespace Ipopt
{
  static const Index dbg_verbosity = 0;

  IpoptAlgorithm::IpoptAlgorithm(const SmartPtr<PDSystemSolver>& pd_solver,
                                 const SmartPtr<LineSearch>& line_search,
                                 const SmartPtr<MuUpdate>& mu_update,
                                 const SmartPtr<ConvergenceCheck>& conv_check,
                                 const SmartPtr<IterateInitializer>& iterate_initializer,
                                 const SmartPtr<IterationOutput>& iter_output)
      :
      pd_solver_(pd_solver),
      line_search_(line_search),
      mu_update_(mu_update),
      conv_check_(conv_check),
      iterate_initializer_(iterate_initializer),
      iter_output_(iter_output)
  {
    DBG_ASSERT(IsValid(pd_solver_));
    DBG_ASSERT(IsValid(line_search_));
    DBG_ASSERT(IsValid(mu_update_));
    DBG_ASSERT(IsValid(conv_check_));
    DBG_ASSERT(IsValid(iterate_initializer_));
    DBG_ASSERT(IsValid(iter_output_));
  }

  IpoptAlgorithm::~IpoptAlgorithm()
  {}

  bool IpoptAlgorithm::InitializeImpl(const OptionsList& options,
                                      const std::string& prefix)
  {

    // Read the IpoptAlgorithm options
    // Initialize the Data object
    bool retvalue = IpData().Initialize(Jnlst(),
                                        options, prefix);
    ASSERT_EXCEPTION(retvalue, AlgorithmStrategyObject::FAILED_INITIALIZATION,
                     "the IpIpoptData object failed to initialize.");

    // Initialize the CQ object
    retvalue = IpCq().Initialize(Jnlst(),
                                 options, prefix);
    ASSERT_EXCEPTION(retvalue, AlgorithmStrategyObject::FAILED_INITIALIZATION,
                     "the IpIpoptCalculatedQuantities object failed to initialize.");

    // Initialize the CQ object
    retvalue = IpNLP().Initialize(Jnlst(),
                                  options, prefix);
    ASSERT_EXCEPTION(retvalue, AlgorithmStrategyObject::FAILED_INITIALIZATION,
                     "the IpIpoptNLP object failed to initialize.");

    // Initialize all the strategies
    retvalue = iterate_initializer_->Initialize(Jnlst(), IpNLP(), IpData(),
               IpCq(), options, prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "the iterate_initializer strategy failed to initialize.");

    retvalue = mu_update_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                      options, prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "the mu_update strategy failed to initialize.");

    retvalue = pd_solver_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                      options, prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "the pd_solver strategy failed to initialize.");

    retvalue = line_search_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                        options,prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "the line_search strategy failed to initialize.");

    retvalue = conv_check_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                       options, prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "the conv_check strategy failed to initialize.");

    retvalue = iter_output_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                        options, prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "the iter_output strategy failed to initialize.");

    if (prefix=="resto.") {
      skip_print_problem_stats_ = true;
    }
    else {
      skip_print_problem_stats_ = false;
    }

    return true;
  }

  IpoptAlgorithm::SolverReturn IpoptAlgorithm::Optimize()
  {
    DBG_START_METH("IpoptAlgorithm::Optimize", dbg_verbosity);

    // Initialize the iterates
    InitializeIterates();

    if (!skip_print_problem_stats_) {
      PrintProblemStatistics();
    }

    ConvergenceCheck::ConvergenceStatus conv_status
    = conv_check_->CheckConvergence();

    // main loop
    while (conv_status == ConvergenceCheck::CONTINUE) {
      // Set the Hessian Matrix
      ActualizeHessian();

      // do all the output for this iteration
      OutputIteration();
      IpData().ResetInfo();

      // update the barrier parameter if necessary
      UpdateBarrierParameter();

      // solve the primal-dual system to get the full step
      ComputeSearchDirection();

      // Compute the new iterate
      ComputeAcceptableTrialPoint();
      //ApplyFractionToBoundary();

      // Accept the new iterate
      AcceptTrialPoint();

      IpData().Set_iter_count(IpData().iter_count()+1);

      conv_status  = conv_check_->CheckConvergence();
    }

    OutputIteration();

    if (conv_status == ConvergenceCheck::CONVERGED) {
      return SUCCESS;
    }
    else if (conv_status == ConvergenceCheck::MAXITER_EXCEEDED) {
      return MAXITER_EXCEEDED;
    }

    return FAILED;
  }

  void IpoptAlgorithm::ActualizeHessian()
  {
    // At this point, just compute the exact Hessian
    IpData().Set_W(IpCq().curr_exact_hessian());
  }


  void IpoptAlgorithm::UpdateBarrierParameter()
  {
    Jnlst().Printf(J_DETAILED, J_MAIN, "\n**************************************************\n");
    Jnlst().Printf(J_DETAILED, J_MAIN, "*** Update Barrier Parameter for Iteration %d:", IpData().iter_count());
    Jnlst().Printf(J_DETAILED, J_MAIN, "\n**************************************************\n\n");
    mu_update_->UpdateBarrierParameter();

    Jnlst().Printf(J_DETAILED, J_MAIN, "Barrier Parameter: %e\n", IpData().curr_mu());

  }

  void IpoptAlgorithm::ComputeSearchDirection()
  {
    DBG_START_METH("IpoptAlgorithm::ComputeSearchDirection", dbg_verbosity);

    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "\n**************************************************\n");
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "*** Solving the Primal Dual System for Iteration %d:",
                   IpData().iter_count());
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "\n**************************************************\n\n");

    if (IpData().HaveDeltas()) {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "No need to compute search direction - it has already been computed.\n");
      return;
    }

    SmartPtr<const Vector> rhs_grad_lag_x  = IpCq().curr_grad_lag_x();
    SmartPtr<const Vector> rhs_grad_lag_s  = IpCq().curr_grad_lag_s();
    SmartPtr<const Vector> rhs_c = IpCq().curr_c();
    SmartPtr<const Vector> rhs_d_minus_s = IpCq().curr_d_minus_s();
    SmartPtr<const Vector> rhs_rel_compl_x_L = IpCq().curr_relaxed_compl_x_L();
    SmartPtr<const Vector> rhs_rel_compl_x_U = IpCq().curr_relaxed_compl_x_U();
    SmartPtr<const Vector> rhs_rel_compl_s_L = IpCq().curr_relaxed_compl_s_L();
    SmartPtr<const Vector> rhs_rel_compl_s_U = IpCq().curr_relaxed_compl_s_U();

    DBG_PRINT_VECTOR(2, "rhs_grad_lag_x", *rhs_grad_lag_x);

    pd_solver_->Solve(-1.0, 0.0,
                      *rhs_grad_lag_x,
                      *rhs_grad_lag_s,
                      *rhs_c,
                      *rhs_d_minus_s,
                      *rhs_rel_compl_x_L,
                      *rhs_rel_compl_x_U,
                      *rhs_rel_compl_s_L,
                      *rhs_rel_compl_s_U,
                      *IpData().NonConst_delta_x(),
                      *IpData().NonConst_delta_s(),
                      *IpData().NonConst_delta_y_c(),
                      *IpData().NonConst_delta_y_d(),
                      *IpData().NonConst_delta_z_L(),
                      *IpData().NonConst_delta_z_U(),
                      *IpData().NonConst_delta_v_L(),
                      *IpData().NonConst_delta_v_U()
                     );

    Jnlst().Printf(J_MOREVECTOR, J_MAIN,
                   "*** Step Calculated for Iteration: %d\n",
                   IpData().iter_count());
    Jnlst().PrintVector(J_MOREVECTOR, J_MAIN, "delta_x", *IpData().delta_x());
    Jnlst().PrintVector(J_MOREVECTOR, J_MAIN, "delta_s", *IpData().delta_s());
    Jnlst().PrintVector(J_MOREVECTOR, J_MAIN, "delta_y_c", *IpData().delta_y_c());
    Jnlst().PrintVector(J_MOREVECTOR, J_MAIN, "delta_y_d", *IpData().delta_y_d());
    Jnlst().PrintVector(J_MOREVECTOR, J_MAIN, "delta_z_L", *IpData().delta_z_L());
    Jnlst().PrintVector(J_MOREVECTOR, J_MAIN, "delta_z_U", *IpData().delta_z_U());
    Jnlst().PrintVector(J_MOREVECTOR, J_MAIN, "delta_v_L", *IpData().delta_v_L());
    Jnlst().PrintVector(J_MOREVECTOR, J_MAIN, "delta_v_U", *IpData().delta_v_U());
  }

  void IpoptAlgorithm::ComputeAcceptableTrialPoint()
  {
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "\n**************************************************\n");
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "*** Finding Acceptable Trial Point for Iteration %d:",
                   IpData().iter_count());
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "\n**************************************************\n\n");
    line_search_->FindAcceptableTrialPoint();
  }

  void IpoptAlgorithm::OutputIteration()
  {
    iter_output_->WriteOutput();
  }

  void IpoptAlgorithm::InitializeIterates()
  {
    DBG_START_METH("IpoptAlgorithm::InitializeIterates", 0);

    iterate_initializer_->SetInitialIterates();
  }

  void IpoptAlgorithm::AcceptTrialPoint()
  {
    // Adjust the bounds if necessary
    Index adjusted_slacks = IpCq().AdjustedTrialSlacks();
    if (adjusted_slacks>0) {
      IpCq().ResetAdjustedTrialSlacks();
      if (adjusted_slacks==1) {
        Jnlst().Printf(J_SUMMARY, J_MAIN,
                       "%d Slack too small, adjusting variable bound\n",
                       adjusted_slacks);
      }
      else {
        Jnlst().Printf(J_SUMMARY, J_MAIN,
                       "%d Slacks too small, adjusting variable bounds\n",
                       adjusted_slacks);
      }
      if (Jnlst().ProduceOutput(J_VECTOR, J_MAIN)) {
        Jnlst().PrintVector(J_VECTOR, J_MAIN, "old_x_L", *IpNLP().x_L());
        Jnlst().PrintVector(J_VECTOR, J_MAIN, "old_x_U", *IpNLP().x_U());
        Jnlst().PrintVector(J_VECTOR, J_MAIN, "old_d_L", *IpNLP().d_L());
        Jnlst().PrintVector(J_VECTOR, J_MAIN, "old_d_U", *IpNLP().d_U());
      }

      SmartPtr<Vector> new_x_l = IpNLP().x_L()->MakeNew();
      IpNLP().Px_L()->TransMultVector(1.0, *IpData().trial_x(),
                                      0.0, *new_x_l);
      new_x_l->Axpy(-1.0, *IpCq().trial_slack_x_L());

      SmartPtr<Vector> new_x_u = IpNLP().x_U()->MakeNew();
      IpNLP().Px_U()->TransMultVector(1.0, *IpData().trial_x(),
                                      0.0, *new_x_u);
      new_x_u->Axpy(1.0, *IpCq().trial_slack_x_U());

      SmartPtr<Vector> new_d_l = IpNLP().d_L()->MakeNew();
      IpNLP().Pd_L()->TransMultVector(1.0, *IpData().trial_s(),
                                      0.0, *new_d_l);
      new_d_l->Axpy(-1.0, *IpCq().trial_slack_s_L());

      SmartPtr<Vector> new_d_u = IpNLP().d_U()->MakeNew();
      IpNLP().Pd_U()->TransMultVector(1.0, *IpData().trial_s(),
                                      0.0, *new_d_u);
      new_d_u->Axpy(1.0, *IpCq().trial_slack_s_U());

      IpNLP().AdjustVariableBounds(*new_x_l, *new_x_u, *new_d_l, *new_d_u);

      if (Jnlst().ProduceOutput(J_VECTOR, J_MAIN)) {
        Jnlst().PrintVector(J_VECTOR, J_MAIN, "new_x_L", *IpNLP().x_L());
        Jnlst().PrintVector(J_VECTOR, J_MAIN, "new_x_U", *IpNLP().x_U());
        Jnlst().PrintVector(J_VECTOR, J_MAIN, "new_d_L", *IpNLP().d_L());
        Jnlst().PrintVector(J_VECTOR, J_MAIN, "new_d_U", *IpNLP().d_U());
      }
    }

    // Accept the step
    IpData().AcceptTrialPoint();
  }

  void IpoptAlgorithm::PrintProblemStatistics()
  {
    if (!Jnlst().ProduceOutput(J_SUMMARY, J_STATISTICS)) {
      // nothing to print
      return;
    }

    SmartPtr<const Vector> x = IpData().curr_x();
    SmartPtr<const Vector> x_L = IpNLP().x_L();
    SmartPtr<const Vector> x_U = IpNLP().x_U();
    SmartPtr<const Matrix> Px_L = IpNLP().Px_L();
    SmartPtr<const Matrix> Px_U = IpNLP().Px_U();

    Index nx_tot, nx_only_lower, nx_both, nx_only_upper;
    calc_number_of_bounds(*IpData().curr_x(), *IpNLP().x_L(), *IpNLP().x_U(),
                          *IpNLP().Px_L(), *IpNLP().Px_U(),
                          nx_tot, nx_only_lower, nx_both, nx_only_upper);

    Index ns_tot, ns_only_lower, ns_both, ns_only_upper;
    calc_number_of_bounds(*IpData().curr_s(), *IpNLP().d_L(), *IpNLP().d_U(),
                          *IpNLP().Pd_L(), *IpNLP().Pd_U(),
                          ns_tot, ns_only_lower, ns_both, ns_only_upper);

    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "Total number of variables............................: %8d\n",nx_tot);
    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "                     variables with only lower bounds: %8d\n",
                   nx_only_lower);
    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "                variables with lower and upper bounds: %8d\n",nx_both);
    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "                     variables with only upper bounds: %8d\n",
                   nx_only_upper);
    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "Total number of equality constraints.................: %8d\n",
                   IpData().curr_y_c()->Dim());
    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "Total number of inequality constraints...............: %8d\n",ns_tot);
    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "        inequality constraints with only lower bounds: %8d\n",
                   ns_only_lower);
    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "   inequality constraints with lower and upper bounds: %8d\n",ns_both);
    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "        inequality constraints with only upper bounds: %8d\n\n",
                   ns_only_upper);
  }

  void IpoptAlgorithm::calc_number_of_bounds(
    const Vector& x,
    const Vector& x_L,
    const Vector& x_U,
    const Matrix& Px_L,
    const Matrix& Px_U,
    Index& n_tot,
    Index& n_only_lower,
    Index& n_both,
    Index& n_only_upper)
  {
    DBG_START_METH("IpoptAlgorithm::calc_number_of_bounds",
                   dbg_verbosity);

    n_tot = x.Dim();

    SmartPtr<Vector> tmpx = x.MakeNew();
    SmartPtr<Vector> tmpxL = x_L.MakeNew();
    SmartPtr<Vector> tmpxU = x_U.MakeNew();

    tmpxL->Set(-1.);
    tmpxU->Set(2.);
    Px_L.MultVector(1.0, *tmpxL, 0.0, *tmpx);
    Px_U.MultVector(1.0, *tmpxU, 1.0, *tmpx);
    // Now, x has elements
    //  -1 : if component has only lower bound
    //   0 : if component has no bound
    //   1 : if component has both lower and upper bound
    //   2 : if component has only upper bound
    DBG_PRINT_VECTOR(2, "x-indicator", *tmpx);

    SmartPtr<Vector> tmpx0 = x.MakeNew();
    tmpx0->Set(0.);

    SmartPtr<Vector> tmpx2 = x.MakeNew();
    tmpx2->Set(-1.0);
    tmpx2->Axpy(1.0, *tmpx);
    tmpx2->ElementWiseMax(*tmpx0); // tmpx2 is now 1 in those
    // components with only upper bounds
    n_only_upper = (Index)tmpx2->Asum();

    tmpx->Axpy(-2., *tmpx2);       // now make all those entries for
    // only upper bounds zero in tmpx

    tmpx2->Copy(*tmpx);
    tmpx2->ElementWiseMax(*tmpx0); // tmpx2 is now 1 in those
    // components with both bounds
    n_both = (Index)tmpx2->Asum();

    tmpx->Axpy(-1., *tmpx2);
    tmpx->ElementWiseMin(*tmpx);   // tmpx is now -1 in those with only
    // lower bounds
    n_only_lower = (Index)tmpx->Asum();

  }

} // namespace Ipopt
