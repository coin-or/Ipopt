// Copyright (C) 2004, International BusinDess Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpPDFullSpaceSolver.hpp"
#include "IpDebug.hpp"

namespace Ipopt
{

  DBG_SET_VERBOSITY(0);

  PDFullSpaceSolver::PDFullSpaceSolver(AugSystemSolver& augSysSolver)
      :
      PDSystemSolver(),
      augSysSolver_(&augSysSolver),
      dummy_cache_(1)
  {
    DBG_START_METH("PDFullSpaceSolver::PDFullSpaceSolver(SmartPtr<AugSystemSolver> augSysSolver)",dbg_verbosity);
  }

  PDFullSpaceSolver::~PDFullSpaceSolver()
  {
    DBG_START_METH("PDFullSpaceSolver::~PDFullSpaceSolver()",dbg_verbosity);
  }


  bool PDFullSpaceSolver::InitializeImpl(const OptionsList& options,
                                         const std::string& prefix)
  {
    Index ivalue;
    Number value;

    // Check for the algorithm options
    if (options.GetIntegerValue("num_min_iter_ref", ivalue, prefix)) {
      ASSERT_EXCEPTION(ivalue >= 0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"num_min_iter_ref\": This value must be larger than or equal to 0");
      num_min_iter_ref_ = ivalue;
    }
    else {
      num_min_iter_ref_ = 1;
    }

    if (options.GetIntegerValue("num_max_iter_ref", ivalue, prefix)) {
      ASSERT_EXCEPTION(ivalue >= num_min_iter_ref_, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"num_max_iter_ref\": This value must be larger than or equal to num_min_iter_ref (default 1)");
      num_max_iter_ref_ = ivalue;
    }
    else {
      num_max_iter_ref_ = 10;
    }

    if (options.GetNumericValue("delta_regu_max", value, prefix)) {
      ASSERT_EXCEPTION(value > 0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"theta_min_fact\": This value must be larger than 0.");
      delta_regu_max_ = value;
    }
    else {
      delta_regu_max_ = 1e40;
    }

    if (options.GetNumericValue("residual_ratio_max", value, prefix)) {
      ASSERT_EXCEPTION(value > 0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"residual_ratio_max\": This value must be larger than 0.");
      residual_ratio_max_ = value;
    }
    else {
      residual_ratio_max_ = 1e-10;
    }

    if (options.GetNumericValue("residual_ratio_singular", value, prefix)) {
      ASSERT_EXCEPTION(value > residual_ratio_max_, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"residual_ratio_singular\": This value must be larger than residual_ratio_max.");
      residual_ratio_singular_ = value;
    }
    else {
      residual_ratio_singular_ = 1e-5;
    }

    if (options.GetNumericValue("residual_improvement_factor", value, prefix)) {
      ASSERT_EXCEPTION(value > 0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"residual_improvement_factor\": This value must be larger than 0.");
      residual_improvement_factor_ = value;
    }
    else {
      residual_improvement_factor_ = 0.999999999;
    }

    // Reset internal flags and data
    augsys_improved_ = false;
    delta_x_curr_=0.;
    delta_s_curr_=0.;
    delta_c_curr_=0.;
    delta_d_curr_=0.;
    delta_x_last_=0.;
    delta_s_last_=0.;
    delta_c_last_=0.;
    delta_d_last_=0.;

    return augSysSolver_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                     options, prefix);
  }

  void PDFullSpaceSolver::Solve(Number alpha,
                                Number beta,
                                const Vector& rhs_x,
                                const Vector& rhs_s,
                                const Vector& rhs_c,
                                const Vector& rhs_d,
                                const Vector& rhs_zL,
                                const Vector& rhs_zU,
                                const Vector& rhs_vL,
                                const Vector& rhs_vU,
                                Vector& res_x,
                                Vector& res_s,
                                Vector& res_c,
                                Vector& res_d,
                                Vector& res_zL,
                                Vector& res_zU,
                                Vector& res_vL,
                                Vector& res_vU,
                                bool allow_inexact)
  {
    DBG_START_METH("PDFullSpaceSolver::Solve",dbg_verbosity);

    DBG_PRINT_VECTOR(2, "rhs_x", rhs_x);
    DBG_PRINT_VECTOR(2, "rhs_s", rhs_s);
    DBG_PRINT_VECTOR(2, "rhs_c", rhs_c);
    DBG_PRINT_VECTOR(2, "rhs_d", rhs_d);
    DBG_PRINT_VECTOR(2, "rhs_zL", rhs_zL);
    DBG_PRINT_VECTOR(2, "rhs_zU", rhs_zU);
    DBG_PRINT_VECTOR(2, "rhs_vL", rhs_vL);
    DBG_PRINT_VECTOR(2, "rhs_vU", rhs_vU);

    // if beta is nonzero, keep a copy of the incoming values in res_ */
    SmartPtr<Vector> copy_res_x;
    SmartPtr<Vector> copy_res_s;
    SmartPtr<Vector> copy_res_c;
    SmartPtr<Vector> copy_res_d;
    SmartPtr<Vector> copy_res_zL;
    SmartPtr<Vector> copy_res_zU;
    SmartPtr<Vector> copy_res_vL;
    SmartPtr<Vector> copy_res_vU;
    if (beta != 0.) {
      copy_res_x = res_x.MakeNew();
      copy_res_s = res_s.MakeNew();
      copy_res_c = res_c.MakeNew();
      copy_res_d = res_d.MakeNew();
      copy_res_zL = res_zL.MakeNew();
      copy_res_zU = res_zU.MakeNew();
      copy_res_vL = res_vL.MakeNew();
      copy_res_vU = res_vU.MakeNew();

      copy_res_x->Copy(res_x);
      copy_res_s->Copy(res_s);
      copy_res_c->Copy(res_c);
      copy_res_d->Copy(res_d);
      copy_res_zL->Copy(res_zL);
      copy_res_zU->Copy(res_zU);
      copy_res_vL->Copy(res_vL);
      copy_res_vU->Copy(res_vU);
    }

    // Receive data about matrix
    SmartPtr<const Vector> x = IpData().curr_x();
    SmartPtr<const Vector> s = IpData().curr_s();
    SmartPtr<const SymMatrix> W = IpData().W();
    SmartPtr<const Matrix> J_c = IpCq().curr_jac_c();
    SmartPtr<const Matrix> J_d = IpCq().curr_jac_d();
    SmartPtr<const Matrix> Px_L = IpNLP().Px_L();
    SmartPtr<const Matrix> Px_U = IpNLP().Px_U();
    SmartPtr<const Matrix> Pd_L = IpNLP().Pd_L();
    SmartPtr<const Matrix> Pd_U = IpNLP().Pd_U();
    SmartPtr<const Vector> z_L = IpData().curr_z_L();
    SmartPtr<const Vector> z_U = IpData().curr_z_U();
    SmartPtr<const Vector> v_L = IpData().curr_v_L();
    SmartPtr<const Vector> v_U = IpData().curr_v_U();
    SmartPtr<const Vector> slack_x_L = IpCq().curr_slack_x_L();
    SmartPtr<const Vector> slack_x_U = IpCq().curr_slack_x_U();
    SmartPtr<const Vector> slack_s_L = IpCq().curr_slack_s_L();
    SmartPtr<const Vector> slack_s_U = IpCq().curr_slack_s_U();
    SmartPtr<const Vector> sigma_x = IpCq().curr_sigma_x();
    SmartPtr<const Vector> sigma_s = IpCq().curr_sigma_s();
    DBG_PRINT_VECTOR(2, "Sigma_x", *sigma_x);
    DBG_PRINT_VECTOR(2, "Sigma_s", *sigma_s);

    bool done = false;
    // The following flag is set to true, if we asked the linear
    // solver successfully, to improve the quality of the solution in
    // the next solve.
    bool resolve_better_quality = false;
    // the following flag is set to true, if iterative refinement
    // failed, and we want to do try if a modified system is able to
    // remedy that problem
    bool pretend_singular = false;
    bool pretend_singular_last_time = false;

    // Beginning of loop for solving the system (including all
    // modifications for the linear system to ensure good solution
    // quality)
    while (!done) {

      bool solve_retval =
        SolveOnce(resolve_better_quality, pretend_singular,
                  *W, *J_c, *J_d, *Px_L, *Px_U, *Pd_L, *Pd_U, *z_L, *z_U,
                  *v_L, *v_U, *slack_x_L, *slack_x_U, *slack_s_L, *slack_s_U,
                  *sigma_x, *sigma_s, 1., 0., rhs_x, rhs_s, rhs_c, rhs_d,
                  rhs_zL, rhs_zU, rhs_vL, rhs_vU, res_x, res_s, res_c, res_d,
                  res_zL, res_zU, res_vL, res_vU);
      resolve_better_quality = false;
      pretend_singular = false;

      if (!solve_retval) {
        // If system seems not to be solvable, we set the search
        // direction to zero, and hope that the line search will take
        // care of this (e.g. call the restoration phase).  ToDo: We
        // might want to use a more explicit cue later.
        res_x.Set(0.);
        res_s.Set(0.);
        res_c.Set(0.);
        res_d.Set(0.);
        res_zL.Set(0.);
        res_zU.Set(0.);
        res_vL.Set(0.);
        res_vU.Set(0.);
        return;
      }

      // Get space for the residual
      SmartPtr<Vector> resid_x = res_x.MakeNew();
      SmartPtr<Vector> resid_s = res_s.MakeNew();
      SmartPtr<Vector> resid_c = res_c.MakeNew();
      SmartPtr<Vector> resid_d = res_d.MakeNew();
      SmartPtr<Vector> resid_zL = res_zL.MakeNew();
      SmartPtr<Vector> resid_zU = res_zU.MakeNew();
      SmartPtr<Vector> resid_vL = res_vL.MakeNew();
      SmartPtr<Vector> resid_vU = res_vU.MakeNew();

      ComputeResiduals(*W, *J_c, *J_d, *Px_L, *Px_U, *Pd_L, *Pd_U,
                       *z_L, *z_U, *v_L, *v_U, *slack_x_L, *slack_x_U,
                       *slack_s_L, *slack_s_U, *sigma_x, *sigma_s,
                       alpha, beta, rhs_x, rhs_s, rhs_c, rhs_d,
                       rhs_zL, rhs_zU, rhs_vL, rhs_vU, res_x, res_s,
                       res_c, res_d,
                       res_zL, res_zU, res_vL, res_vU, *resid_x, *resid_s,
                       *resid_c, *resid_d,
                       *resid_zL, *resid_zU, *resid_vL, *resid_vU);

      Number residual_ratio =
        ComputeResidualRatio(rhs_x, rhs_s, rhs_c, rhs_d,
                             rhs_zL, rhs_zU, rhs_vL, rhs_vU,
                             res_x, res_s, res_c, res_d,
                             res_zL, res_zU, res_vL, res_vU,
                             *resid_x, *resid_s, *resid_c, *resid_d,
                             *resid_zL, *resid_zU, *resid_vL, *resid_vU);
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "residual_ratio = %e\n", residual_ratio);
      Number residual_ratio_old = residual_ratio;

      // Beginning of loop for iterative refinement
      Index num_iter_ref = 0;
      bool quit_refinement = false;
      while (!allow_inexact && !quit_refinement &&
             (num_iter_ref < num_min_iter_ref_ ||
              residual_ratio > residual_ratio_max_) ) {

        // To the next back solve
        solve_retval =
          SolveOnce(resolve_better_quality, false,
                    *W, *J_c, *J_d, *Px_L, *Px_U, *Pd_L, *Pd_U, *z_L, *z_U,
                    *v_L, *v_U, *slack_x_L, *slack_x_U, *slack_s_L, *slack_s_U,
                    *sigma_x, *sigma_s, -1., 1., *resid_x, *resid_s, *resid_c, *resid_d,
                    *resid_zL, *resid_zU, *resid_vL, *resid_vU, res_x, res_s, res_c, res_d,
                    res_zL, res_zU, res_vL, res_vU);
        DBG_ASSERT(solve_retval);

        ComputeResiduals(*W, *J_c, *J_d, *Px_L, *Px_U, *Pd_L, *Pd_U,
                         *z_L, *z_U, *v_L, *v_U, *slack_x_L, *slack_x_U,
                         *slack_s_L, *slack_s_U, *sigma_x, *sigma_s,
                         alpha, beta, rhs_x, rhs_s, rhs_c, rhs_d,
                         rhs_zL, rhs_zU, rhs_vL, rhs_vU, res_x, res_s,
                         res_c, res_d,
                         res_zL, res_zU, res_vL, res_vU, *resid_x, *resid_s,
                         *resid_c, *resid_d,
                         *resid_zL, *resid_zU, *resid_vL, *resid_vU);

        residual_ratio =
          ComputeResidualRatio(rhs_x, rhs_s, rhs_c, rhs_d,
                               rhs_zL, rhs_zU, rhs_vL, rhs_vU,
                               res_x, res_s, res_c, res_d,
                               res_zL, res_zU, res_vL, res_vU,
                               *resid_x, *resid_s, *resid_c, *resid_d,
                               *resid_zL, *resid_zU, *resid_vL, *resid_vU);
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                       "residual_ratio = %e\n", residual_ratio);

        num_iter_ref++;
        // Check if we have to give up on iterative refinement
        if (num_iter_ref>num_min_iter_ref_ &&
            (num_iter_ref>num_max_iter_ref_ ||
             residual_ratio>residual_improvement_factor_*residual_ratio_old)) {

          Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                         "Iterative refinement failed with residual_ratio = %e\n", residual_ratio);
          quit_refinement = true;

          // Pretend singularity only once - if it didn't help, we
          // have to live with what we got so far
          resolve_better_quality = false;
	  DBG_PRINT((1, "pretend_singular = %d\n", pretend_singular));
          if (!pretend_singular_last_time) {
            // First try if we can ask the augmented system solver to
            // improve the quality of the solution (only if that hasn't
            // been done before for this linear system)
            if (!augsys_improved_) {
              Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                             "Asking augmented system solver to improve quality of its solutions.\n");
              augsys_improved_ = augSysSolver_->IncreaseQuality();
              if (augsys_improved_) {
                IpData().Append_info_string("q");
		resolve_better_quality = true;
              }
              else {
                // solver said it cannot improve quality, so let
                // possibly conclude that the current modification is
                // singular
                pretend_singular = true;
              }
            }
            else {
              // we had already asked the solver before to improve the
              // quality of the solution, so let's now pretend that the
              // modification is possibly singular
              pretend_singular = true;
            }
	    pretend_singular_last_time = pretend_singular;
            if (pretend_singular) {
              // let's only conclude that the current linear system
              // including modifications is singular, if the residual is
              // quite bad
              if (residual_ratio < residual_ratio_singular_) {
                pretend_singular = false;
                IpData().Append_info_string("S");
                Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                               "Just accept current solution.\n");
              }
              else {
                IpData().Append_info_string("s");
                Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                               "Pretend that the current system (including modifications) is singular.\n");
              }
            }
          }
          else {
            pretend_singular = false;
	    DBG_PRINT((1,"Resetting pretend_singular to false.\n"));
          }
        }

        residual_ratio_old = residual_ratio;
      } // End of loop for iterative refinement

      done = !(resolve_better_quality) && !(pretend_singular);

    } // End of loop for solving the linear system (incl. modifications)

    // Now that the system has been solved, remember the current
    // perturbation for the next time
    if (delta_x_curr_!=0.) {
      delta_x_last_ = delta_x_curr_;
    }
    if (delta_s_curr_!=0.) {
      delta_s_last_ = delta_s_curr_;
    }
    if (delta_c_curr_!=0.) {
      delta_c_last_ = delta_c_curr_;
      IpData().Append_info_string("L");
    }
    if (delta_s_curr_!=0.) {
      delta_d_last_ = delta_d_curr_;
    }

    // Finally let's assemble the res result vectors
    if (alpha != 0.) {
      res_x.Scal(alpha);
      res_s.Scal(alpha);
      res_c.Scal(alpha);
      res_d.Scal(alpha);
      res_zL.Scal(alpha);
      res_zU.Scal(alpha);
      res_vL.Scal(alpha);
      res_vU.Scal(alpha);
    }

    if (beta != 0.) {
      res_x.Axpy(beta, *copy_res_x);
      res_s.Axpy(beta, *copy_res_s);
      res_c.Axpy(beta, *copy_res_c);
      res_d.Axpy(beta, *copy_res_d);
      res_zL.Axpy(beta, *copy_res_zL);
      res_zU.Axpy(beta, *copy_res_zU);
      res_vL.Axpy(beta, *copy_res_vL);
      res_vU.Axpy(beta, *copy_res_vU);
    }

    // Set some information for iteration summary output
    IpData().Set_info_regu_x(delta_x_curr_);

  }

  bool PDFullSpaceSolver::SolveOnce(bool resolve_unmodified,
                                    bool pretend_singular,
                                    const SymMatrix& W,
                                    const Matrix& J_c,
                                    const Matrix& J_d,
                                    const Matrix& Px_L,
                                    const Matrix& Px_U,
                                    const Matrix& Pd_L,
                                    const Matrix& Pd_U,
                                    const Vector& z_L,
                                    const Vector& z_U,
                                    const Vector& v_L,
                                    const Vector& v_U,
                                    const Vector& slack_x_L,
                                    const Vector& slack_x_U,
                                    const Vector& slack_s_L,
                                    const Vector& slack_s_U,
                                    const Vector& sigma_x,
                                    const Vector& sigma_s,
                                    Number alpha,
                                    Number beta,
                                    const Vector& rhs_x,
                                    const Vector& rhs_s,
                                    const Vector& rhs_c,
                                    const Vector& rhs_d,
                                    const Vector& rhs_zL,
                                    const Vector& rhs_zU,
                                    const Vector& rhs_vL,
                                    const Vector& rhs_vU,
                                    Vector& res_x,
                                    Vector& res_s,
                                    Vector& res_c,
                                    Vector& res_d,
                                    Vector& res_zL,
                                    Vector& res_zU,
                                    Vector& res_vL,
                                    Vector& res_vU)
  {
    // TO DO LIST:
    //
    // 1. decide for reasonable return codes (e.g. fatal error, too
    //    ill-conditioned...)
    // 2. Make constants parameters that can be set from the outside
    // 3. Get Information out of Ipopt structures
    // 4. add heuristic for structurally singular problems
    // 5. see if it makes sense to distinguish delta_x and delta_s,
    //    or delta_c and delta_d
    // 6. increase pivot tolerance if number of get evals so too small
    DBG_START_METH("PDFullSpaceSolver::SolveOnce",dbg_verbosity);

    // Compute the right hand side for the augmented system formulation
    SmartPtr<Vector> augRhs_x = rhs_x.MakeNew();
    augRhs_x->Copy(rhs_x);
    Px_L.AddMSinvZ(1.0, slack_x_L, rhs_zL, *augRhs_x);
    Px_U.AddMSinvZ(-1.0, slack_x_U, rhs_zU, *augRhs_x);
    /*
    AddPSinvZ(1.0, Px_L, slack_x_L, rhs_zL, 1.0, *augRhs_x);
    AddPSinvZ(-1.0, Px_U, slack_x_U, rhs_zU, 1.0, *augRhs_x);
    */

    SmartPtr<Vector> augRhs_s = rhs_s.MakeNew();
    augRhs_s->Copy(rhs_s);
    Pd_L.AddMSinvZ(1.0, slack_s_L, rhs_vL, *augRhs_s);
    Pd_U.AddMSinvZ(-1.0, slack_s_U, rhs_vU, *augRhs_s);
    /*
    AddPSinvZ(1.0, Pd_L, slack_s_L, rhs_vL, 1.0, *augRhs_s);
    AddPSinvZ(-1.0, Pd_U, slack_s_U, rhs_vU, 1.0, *augRhs_s);
    */

    // Get space into which we can put the solution of the augmented system
    SmartPtr<Vector> sol_x = res_x.MakeNew();
    SmartPtr<Vector> sol_s = res_s.MakeNew();
    SmartPtr<Vector> sol_c = res_c.MakeNew();
    SmartPtr<Vector> sol_d = res_d.MakeNew();

    // Now check whether any data has changed
    std::vector<const TaggedObject*> deps;
    deps.push_back(&W);
    deps.push_back(&J_c);
    deps.push_back(&J_d);
    deps.push_back(&z_L);
    deps.push_back(&z_U);
    deps.push_back(&v_L);
    deps.push_back(&v_U);
    deps.push_back(&slack_x_L);
    deps.push_back(&slack_x_U);
    deps.push_back(&slack_s_L);
    deps.push_back(&slack_s_U);
    deps.push_back(&sigma_x);
    deps.push_back(&sigma_s);
    void* dummy;
    bool uptodate = dummy_cache_.GetCachedResult(dummy, deps);
    if (!uptodate) {
      dummy_cache_.AddCachedResult(dummy, deps);
      augsys_improved_ = false;
    }
    // improve_current_solution can only be true, if that system has
    // been solved before
    DBG_ASSERT((!resolve_unmodified && !pretend_singular) || uptodate);

    ESymSolverStatus retval;
    if (uptodate && !pretend_singular) {
      // No need to go throught the pain of finding the appropriate
      // values for the deltas, because the matrix hasn't changed since
      // the last call.  So, just call the Solve Method
      retval = augSysSolver_->Solve(&W, &sigma_x, delta_x_curr_,
                                    &sigma_s, delta_s_curr_, &J_c, NULL,
                                    delta_c_curr_, &J_d, NULL, delta_d_curr_,
                                    *augRhs_x, *augRhs_s, rhs_c, rhs_d,
                                    *sol_x, *sol_s, *sol_c,
                                    *sol_d, false, 0);
      assert(retval==SYMSOLVER_SUCCESS); //TODO make return code
    }
    else {
      Index numberOfEVals=rhs_c.Dim()+rhs_d.Dim();
      // counter for the number of trial evaluations
      // (ToDo is not at the corrent place)
      Index count = 0;
      // Flag indicating if instead of the first
      // solve we want to pretend that the system is
      // singluar
      if (true || !pretend_singular) {
        // First try, if no modification is necessary to obtain correct inertia
        delta_x_curr_=0.;
        delta_s_curr_=0.;
        delta_c_curr_=0.;
        delta_d_curr_=0.;
      }

      retval = SYMSOLVER_SINGULAR;
      bool fail = false;

      while (retval!= SYMSOLVER_SUCCESS && !fail) {

        if (pretend_singular) {
          retval = SYMSOLVER_SINGULAR;
          pretend_singular = false;
        }
        else {
          count++;
          Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                         "Solving system with delta_x=%e delta_s=%e\n                    delta_c=%e delta_d=%e\n",
                         delta_x_curr_, delta_s_curr_,
                         delta_c_curr_, delta_d_curr_);
          retval = augSysSolver_->Solve(&W, &sigma_x, delta_x_curr_,
                                        &sigma_s, delta_s_curr_, &J_c, NULL,
                                        delta_c_curr_, &J_d, NULL, delta_d_curr_,
                                        *augRhs_x, *augRhs_s, rhs_c, rhs_d,
                                        *sol_x, *sol_s, *sol_c,
                                        *sol_d, true, numberOfEVals);
        }
        assert(retval!=SYMSOLVER_FATAL_ERROR); //TODO make return code
        if (retval==SYMSOLVER_SINGULAR && delta_c_curr_==0.) {
          // If the matrix is singular and delta_c is not yet nonzero,
          // increase both delta_c and delta_d
          delta_c_curr_ = 1e-8; // TODO: Make parameter
          delta_d_curr_ = 1e-8; // TODO Make parameter
        }
        else if (retval==SYMSOLVER_WRONG_INERTIA || retval==SYMSOLVER_SINGULAR) {
          if (delta_x_curr_ == 0) {
            if (delta_x_last_ == 0.) {
              delta_x_curr_ = 1e-4;  //TODO Parameter
            }
            else {
              delta_x_curr_ = Max(1e-20,delta_x_last_/3.);  //TODO Parameter
            }
          }
          else {
            if (delta_x_last_ == 0. || 1e5*delta_x_last_<delta_x_curr_) {
              delta_x_curr_ = 100.*delta_x_curr_;  //TODO Parameter
            }
            else {
              delta_x_curr_ = 8.*delta_x_curr_;  //TODO Parameter
            }
          }
          if (delta_x_curr_ > delta_regu_max_) {
            // Give up trying to solve the linear system
            return false;
          }
          delta_s_curr_ = delta_x_curr_;
        }
      } // while (retval!=SYMSOLVER_SUCCESS && !fail) {

      // Some output
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Number of trial factorizations performed: %d\n",
                     count);
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Perturbation parameters: delta_x=%e delta_s=%e\n                         delta_c=%e delta_d=%e\n",
                     delta_x_curr_,delta_s_curr_,delta_c_curr_,delta_d_curr_);
    }

    // Compute the remaining sol Vectors
    SmartPtr<Vector> sol_zL = res_zL.MakeNew();
    SmartPtr<Vector> sol_zU = res_zU.MakeNew();
    SmartPtr<Vector> sol_vL = res_vL.MakeNew();
    SmartPtr<Vector> sol_vU = res_vU.MakeNew();

    Px_L.SinvBlrmZMTdBr(-1., slack_x_L, rhs_zL, z_L, *sol_x, *sol_zL);
    Px_U.SinvBlrmZMTdBr(1., slack_x_U, rhs_zU, z_U, *sol_x, *sol_zU);
    Pd_L.SinvBlrmZMTdBr(-1., slack_s_L, rhs_vL, v_L, *sol_s, *sol_vL);
    Pd_U.SinvBlrmZMTdBr(1., slack_s_U, rhs_vU, v_U, *sol_s, *sol_vU);

    /*
    SinvBlrmZPTdBr(-1., slack_x_L, rhs_zL, z_L, Px_L, *sol_x, *sol_zL);
    SinvBlrmZPTdBr(1., slack_x_U, rhs_zU, z_U, Px_U, *sol_x, *sol_zU);
    SinvBlrmZPTdBr(-1., slack_s_L, rhs_vL, v_L, Pd_L, *sol_s, *sol_vL);
    SinvBlrmZPTdBr(1., slack_s_U, rhs_vU, v_U, Pd_U, *sol_s, *sol_vU);
    */

    // Finally let's assemble the res result vectors
    AxpBy(alpha, *sol_x, beta, res_x);
    AxpBy(alpha, *sol_s, beta, res_s);
    AxpBy(alpha, *sol_c, beta, res_c);
    AxpBy(alpha, *sol_d, beta, res_d);
    AxpBy(alpha, *sol_zL, beta, res_zL);
    AxpBy(alpha, *sol_zU, beta, res_zU);
    AxpBy(alpha, *sol_vL, beta, res_vL);
    AxpBy(alpha, *sol_vU, beta, res_vU);

    return true;
  }

  void PDFullSpaceSolver::ComputeResiduals(
    const SymMatrix& W,
    const Matrix& J_c,
    const Matrix& J_d,
    const Matrix& Px_L,
    const Matrix& Px_U,
    const Matrix& Pd_L,
    const Matrix& Pd_U,
    const Vector& z_L,
    const Vector& z_U,
    const Vector& v_L,
    const Vector& v_U,
    const Vector& slack_x_L,
    const Vector& slack_x_U,
    const Vector& slack_s_L,
    const Vector& slack_s_U,
    const Vector& sigma_x,
    const Vector& sigma_s,
    Number alpha,
    Number beta,
    const Vector& rhs_x,
    const Vector& rhs_s,
    const Vector& rhs_c,
    const Vector& rhs_d,
    const Vector& rhs_zL,
    const Vector& rhs_zU,
    const Vector& rhs_vL,
    const Vector& rhs_vU,
    const Vector& res_x,
    const Vector& res_s,
    const Vector& res_c,
    const Vector& res_d,
    const Vector& res_zL,
    const Vector& res_zU,
    const Vector& res_vL,
    const Vector& res_vU,
    Vector& resid_x,
    Vector& resid_s,
    Vector& resid_c,
    Vector& resid_d,
    Vector& resid_zL,
    Vector& resid_zU,
    Vector& resid_vL,
    Vector& resid_vU)
  {
    DBG_START_METH("PDFullSpaceSolver::ComputeResiduals", dbg_verbosity);

    SmartPtr<Vector> tmp;

    // x
    W.MultVector(1., res_x, 0., resid_x);
    //    if (delta_x_curr_!=0.) {
    //      resid_x.Axpy(delta_x_curr_, res_x);
    //    }
    J_c.TransMultVector(1., res_c, 1., resid_x);
    J_d.TransMultVector(1., res_d, 1., resid_x);
    Px_L.MultVector(-1., res_zL, 1., resid_x);
    Px_U.MultVector(1., res_zU, 1., resid_x);
    resid_x.AddTwoVectors(delta_x_curr_, res_x, -1., rhs_x, 1.);
    //    resid_x.Axpy(-1., rhs_x);

    // s
    Pd_U.MultVector(1., res_vU, 0., resid_s);
    Pd_L.MultVector(-1., res_vL, 1., resid_s);
    resid_s.AddTwoVectors(-1., res_d, -1., rhs_s, 1.);
    if (delta_s_curr_!=0.) {
      resid_s.Axpy(delta_s_curr_, res_s);
    }
    
    /* DELE
    resid_s.Axpy(-1., res_d);
    if (delta_s_curr_!=0.) {
      resid_s.Axpy(delta_s_curr_, res_s);
    }
    resid_s.Axpy(-1., rhs_s);
    */

    // c
    J_c.MultVector(1., res_x, 0., resid_c);
    resid_c.AddTwoVectors(-delta_c_curr_, res_c, -1., rhs_c, 1.);
    /* DELE
    if (delta_c_curr_) {
      resid_c.Axpy(-delta_c_curr_, res_c);
    }
    resid_c.Axpy(-1., rhs_c);
    */

    // d
    J_d.MultVector(1., res_x, 0., resid_d);
    resid_d.AddTwoVectors(-1., res_s, -1., rhs_d, 1.);
    if (delta_d_curr_) {
      resid_d.Axpy(-delta_d_curr_, res_d);
    }

    /* DELE
    resid_d.Axpy(-1., res_s);
    if (delta_d_curr_) {
      resid_d.Axpy(-delta_d_curr_, res_d);
    }
    resid_d.Axpy(-1., rhs_d);
    */

    // zL
    resid_zL.Copy(res_zL);
    resid_zL.ElementWiseMultiply(slack_x_L);
    tmp = z_L.MakeNew();
    Px_L.TransMultVector(1., res_x, 0., *tmp);
    tmp->ElementWiseMultiply(z_L);
    resid_zL.AddTwoVectors(1., *tmp, -1., rhs_zL, 1.);
    /* DELE
    resid_zL.Axpy(1., *tmp);
    resid_zL.Axpy(-1., rhs_zL);
    */

    // zU
    resid_zU.Copy(res_zU);
    resid_zU.ElementWiseMultiply(slack_x_U);
    tmp = z_U.MakeNew();
    Px_U.TransMultVector(1., res_x, 0., *tmp);
    tmp->ElementWiseMultiply(z_U);
    resid_zU.AddTwoVectors(-1., *tmp, -1., rhs_zU, 1.);
    /* DELE
    resid_zU.Axpy(-1., *tmp);
    resid_zU.Axpy(-1., rhs_zU);
    */

    // vL
    resid_vL.Copy(res_vL);
    resid_vL.ElementWiseMultiply(slack_s_L);
    tmp = v_L.MakeNew();
    Pd_L.TransMultVector(1., res_s, 0., *tmp);
    tmp->ElementWiseMultiply(v_L);
    resid_vL.AddTwoVectors(1., *tmp, -1., rhs_vL, 1.);
    /* DELE
    resid_vL.Axpy(1., *tmp);
    resid_vL.Axpy(-1., rhs_vL);
    */

    // vU
    resid_vU.Copy(res_vU);
    resid_vU.ElementWiseMultiply(slack_s_U);
    tmp = v_U.MakeNew();
    Pd_U.TransMultVector(1., res_s, 0., *tmp);
    tmp->ElementWiseMultiply(v_U);
    resid_vU.AddTwoVectors(-1., *tmp, -1., rhs_vU, 1.);
    /* DELE
    resid_vU.Axpy(-1., *tmp);
    resid_vU.Axpy(-1., rhs_vU);
    */

    if (Jnlst().ProduceOutput(J_MOREVECTOR, J_LINEAR_ALGEBRA)) {
      Jnlst().PrintVector(J_MOREVECTOR, J_LINEAR_ALGEBRA,"resid_x ", resid_x);
      Jnlst().PrintVector(J_MOREVECTOR, J_LINEAR_ALGEBRA,"resid_s ", resid_s);
      Jnlst().PrintVector(J_MOREVECTOR, J_LINEAR_ALGEBRA,"resid_c ", resid_c);
      Jnlst().PrintVector(J_MOREVECTOR, J_LINEAR_ALGEBRA,"resid_d ", resid_d);
      Jnlst().PrintVector(J_MOREVECTOR, J_LINEAR_ALGEBRA,"resid_zL", resid_zL);
      Jnlst().PrintVector(J_MOREVECTOR, J_LINEAR_ALGEBRA,"resid_zU", resid_zU);
      Jnlst().PrintVector(J_MOREVECTOR, J_LINEAR_ALGEBRA,"resid_vL", resid_vL);
      Jnlst().PrintVector(J_MOREVECTOR, J_LINEAR_ALGEBRA,"resid_vU", resid_vU);
    }

    if (Jnlst().ProduceOutput(J_MOREDETAILED, J_LINEAR_ALGEBRA)) {
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "max-norm resid_x  %e\n", resid_x.Amax());
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "max-norm resid_s  %e\n", resid_s.Amax());
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "max-norm resid_c  %e\n", resid_c.Amax());
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "max-norm resid_d  %e\n", resid_d.Amax());
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "max-norm resid_zL %e\n", resid_zL.Amax());
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "max-norm resid_zU %e\n", resid_zU.Amax());
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "max-norm resid_vL %e\n", resid_vL.Amax());
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "max-norm resid_vU %e\n", resid_vU.Amax());
    }
  }

  Number PDFullSpaceSolver::ComputeResidualRatio(
    const Vector& rhs_x,
    const Vector& rhs_s,
    const Vector& rhs_c,
    const Vector& rhs_d,
    const Vector& rhs_zL,
    const Vector& rhs_zU,
    const Vector& rhs_vL,
    const Vector& rhs_vU,
    const Vector& res_x,
    const Vector& res_s,
    const Vector& res_c,
    const Vector& res_d,
    const Vector& res_zL,
    const Vector& res_zU,
    const Vector& res_vL,
    const Vector& res_vU,
    const Vector& resid_x,
    const Vector& resid_s,
    const Vector& resid_c,
    const Vector& resid_d,
    const Vector& resid_zL,
    const Vector& resid_zU,
    const Vector& resid_vL,
    const Vector& resid_vU)
  {
    DBG_START_METH("PDFullSpaceSolver::ComputeResidualRatio", dbg_verbosity);

    Number nrm_rhs = rhs_x.Amax();
    nrm_rhs = Max(nrm_rhs, rhs_s.Amax());
    nrm_rhs = Max(nrm_rhs, rhs_c.Amax());
    nrm_rhs = Max(nrm_rhs, rhs_d.Amax());
    nrm_rhs = Max(nrm_rhs, rhs_zL.Amax());
    nrm_rhs = Max(nrm_rhs, rhs_zU.Amax());
    nrm_rhs = Max(nrm_rhs, rhs_vL.Amax());
    nrm_rhs = Max(nrm_rhs, rhs_vU.Amax());

    Number nrm_res = res_x.Amax();
    nrm_res = Max(nrm_res, res_s.Amax());
    nrm_res = Max(nrm_res, res_c.Amax());
    nrm_res = Max(nrm_res, res_d.Amax());
    nrm_res = Max(nrm_res, res_zL.Amax());
    nrm_res = Max(nrm_res, res_zU.Amax());
    nrm_res = Max(nrm_res, res_vL.Amax());
    nrm_res = Max(nrm_res, res_vU.Amax());

    Number nrm_resid = resid_x.Amax();
    nrm_resid = Max(nrm_resid, resid_s.Amax());
    nrm_resid = Max(nrm_resid, resid_c.Amax());
    nrm_resid = Max(nrm_resid, resid_d.Amax());
    nrm_resid = Max(nrm_resid, resid_zL.Amax());
    nrm_resid = Max(nrm_resid, resid_zU.Amax());
    nrm_resid = Max(nrm_resid, resid_vL.Amax());
    nrm_resid = Max(nrm_resid, resid_vU.Amax());

    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "nrm_rhs = %8.2e nrm_sol = %8.2e nrm_resid = %8.2e\n",
                   nrm_rhs, nrm_res, nrm_resid);

    if (nrm_rhs+nrm_res == 0.) {
      return nrm_resid;  // this should be zero
    }
    else {
      // ToDo: determine how to include norm of matrix, and what
      // safeguard to use against incredibly large solution vectors
      Number max_cond = 1e6;
      return nrm_resid/(Min(nrm_res, max_cond*nrm_rhs)+nrm_rhs);
    }
  }

  /*
  void PDFullSpaceSolver::AddPSinvZ(Number alpha, const Matrix& P,
                                    const Vector& S, const Vector& Z,
                                    Number beta, Vector& X)
  {
    SmartPtr<Vector> tmp = S.MakeNew();
    tmp->Copy(Z);
    tmp->ElementWiseDivide(S);
    P.MultVector(alpha, *tmp, beta, X);
  }

  void PDFullSpaceSolver::SinvBlrmZPTdBr(Number alpha, const Vector& S,
                                         const Vector& R, const Vector& Z,
                                         const Matrix& P, const Vector& D,
                                         Vector& X)
  {
    P.TransMultVector(alpha, D, 0., X);
    X.ElementWiseMultiply(Z);
    X.Axpy(1., R);
    X.ElementWiseDivide(S);
  }
  */

  void PDFullSpaceSolver::AxpBy(Number alpha, const Vector& X,
                                Number beta, Vector& Y)
  {
    if (beta==0.) {
      Y.Set(0.);
    }
    else if (beta!=1.) {
      Y.Scal(beta);
    }
    if (alpha!=0.) {
      Y.Axpy(alpha, X);
    }
  }

} // namespace Ipopt
