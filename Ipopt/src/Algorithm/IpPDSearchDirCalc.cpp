// Copyright (C) 2005, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter                 IBM    2005-10-13
//               derived from IpIpoptAlg.cpp

#include "IpPDSearchDirCalc.hpp"

namespace Ipopt
{

#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  PDSearchDirCalculator::PDSearchDirCalculator(const SmartPtr<PDSystemSolver>& pd_solver)
      :
      pd_solver_(pd_solver)
  {
    DBG_START_FUN("PDSearchDirCalculator::PDSearchDirCalculator",
                  dbg_verbosity);
    DBG_ASSERT(IsValid(pd_solver_));
  }

  PDSearchDirCalculator::~PDSearchDirCalculator()
  {
    DBG_START_FUN("PDSearchDirCalculator::~PDSearchDirCalculator()",
                  dbg_verbosity);
  }

  void PDSearchDirCalculator::RegisterOptions(const SmartPtr<RegisteredOptions>& roptions)
  {
    roptions->SetRegisteringCategory("Step Calculation");
    roptions->AddStringOption2(
      "fast_step_computation",
      "Indicates if the linear system should be solved quickly.",
      "no",
      "no", "Verify solution of linear system by computing residuals.",
      "yes", "Trust that linear systems are solved well.",
      "If set to yes, the algorithm assumes that the linear system that is "
      "solved to obtain the search direction, is solved sufficiently well. "
      "In that case, no residuals are computed, and the computation of the "
      "search direction is a little faster.");
  }


  bool PDSearchDirCalculator::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    options.GetBoolValue("fast_step_computation", fast_step_computation_,
                         prefix);
    options.GetBoolValue("mehrotra_algorithm", mehrotra_algorithm_, prefix);
    return pd_solver_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                  options, prefix);
  }

  bool PDSearchDirCalculator::ComputeSearchDirection()
  {
    DBG_START_METH("PDSearchDirCalculator::ComputeSearchDirection",
                   dbg_verbosity);

    bool improve_solution = false;
    if (IpData().HaveDeltas()) {
      improve_solution = true;
    }

    bool retval;
    if (improve_solution && fast_step_computation_) {
      retval = true;
    }
    else {
      SmartPtr<IteratesVector> rhs = IpData().curr()->MakeNewContainer();
      rhs->Set_x(*IpCq().curr_grad_lag_with_damping_x());
      rhs->Set_s(*IpCq().curr_grad_lag_with_damping_s());
      rhs->Set_y_c(*IpCq().curr_c());
      rhs->Set_y_d(*IpCq().curr_d_minus_s());
      Index nbounds = IpNLP().x_L()->Dim()+ IpNLP().x_U()->Dim() +
                      IpNLP().d_L()->Dim()+ IpNLP().d_U()->Dim();
      if (nbounds>0 && mehrotra_algorithm_) {
        // set up the right hand side a la Mehrotra
        DBG_ASSERT(IpData().HaveAffineDeltas());
        DBG_ASSERT(!IpData().HaveDeltas());
        const SmartPtr<const IteratesVector> delta_aff = IpData().delta_aff();

        SmartPtr<Vector> tmpvec = delta_aff->z_L()->MakeNew();
        IpNLP().Px_L()->TransMultVector(1., *delta_aff->x(), 0., *tmpvec);
        tmpvec->ElementWiseMultiply(*delta_aff->z_L());
        tmpvec->Axpy(1., *IpCq().curr_relaxed_compl_x_L());
        rhs->Set_z_L(*tmpvec);

        tmpvec = delta_aff->z_U()->MakeNew();
        IpNLP().Px_U()->TransMultVector(-1., *delta_aff->x(), 0., *tmpvec);
        tmpvec->ElementWiseMultiply(*delta_aff->z_U());
        tmpvec->Axpy(1., *IpCq().curr_relaxed_compl_x_U());
        rhs->Set_z_U(*tmpvec);

        tmpvec = delta_aff->v_L()->MakeNew();
        IpNLP().Pd_L()->TransMultVector(1., *delta_aff->s(), 0., *tmpvec);
        tmpvec->ElementWiseMultiply(*delta_aff->v_L());
        tmpvec->Axpy(1., *IpCq().curr_relaxed_compl_s_L());
        rhs->Set_v_L(*tmpvec);

        tmpvec = delta_aff->v_U()->MakeNew();
        IpNLP().Pd_U()->TransMultVector(-1., *delta_aff->s(), 0., *tmpvec);
        tmpvec->ElementWiseMultiply(*delta_aff->v_U());
        tmpvec->Axpy(1., *IpCq().curr_relaxed_compl_s_U());
        rhs->Set_v_U(*tmpvec);
      }
      else {
        rhs->Set_z_L(*IpCq().curr_relaxed_compl_x_L());
        rhs->Set_z_U(*IpCq().curr_relaxed_compl_x_U());
        rhs->Set_v_L(*IpCq().curr_relaxed_compl_s_L());
        rhs->Set_v_U(*IpCq().curr_relaxed_compl_s_U());
      }

      DBG_PRINT_VECTOR(2, "rhs", *rhs);

      // Get space for the search direction
      SmartPtr<IteratesVector> delta =
        IpData().curr()->MakeNewIteratesVector(true);

      if (improve_solution) {
        // We can probably avoid copying and scaling...
        delta->AddOneVector(-1., *IpData().delta(), 0.);
      }

      bool& allow_inexact = fast_step_computation_;
      retval = pd_solver_->Solve(-1.0, 0.0, *rhs, *delta, allow_inexact,
                                 improve_solution);
      if (retval) {
        // Store the search directions in the IpData object
        IpData().set_delta(delta);
      }
    }
    return retval;
  }

} // namespace Ipopt
