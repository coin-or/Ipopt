// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2008-09-19

#include "IpInexactNormalTerminationTester.hpp"
#include "IpBlas.hpp"

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

  InexactNormalTerminationTester::InexactNormalTerminationTester()
  {
    DBG_START_METH("InexactNormalTerminationTester::InexactNormalTerminationTester",dbg_verbosity);
  }

  InexactNormalTerminationTester::~InexactNormalTerminationTester()
  {
    DBG_START_METH("InexactNormalTerminationTester::~InexactNormalTerminationTester()",dbg_verbosity);
  }

  void InexactNormalTerminationTester::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddLowerBoundedNumberOption(
      "inexact_normal_tol",
      "Desired relative residual tolerance for iterative solver during normal step computation.",
      0.0, true, 1e-3,
      "");
    roptions->AddLowerBoundedIntegerOption(
      "inexact_normal_max_iter",
      "Maximal number of iterative solver iterations during normal step computation",
      0, 200,
      "");
  }


  bool InexactNormalTerminationTester::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    options.GetNumericValue("inexact_normal_tol", inexact_normal_tol_, prefix);
    options.GetIntegerValue("inexact_normal_max_iter",
                            inexact_normal_max_iter_, prefix);

    std::string inexact_linear_system_scaling;
    options.GetStringValue("inexact_linear_system_scaling",
                           inexact_linear_system_scaling, prefix);
    if (inexact_linear_system_scaling=="slack-based") {
      requires_scaling_ = true;
    }
    else {
      requires_scaling_ = false;
    }

    c_Avc_norm_cauchy_ = -1;
    return true;
  }

  bool InexactNormalTerminationTester::InitializeSolve()
  {
    DBG_START_METH("InexactNormalTerminationTester::InitializeSolve",
                   dbg_verbosity);


    return true;
  }

  InexactNormalTerminationTester::ETerminationTest
  InexactNormalTerminationTester::
  TestTermination(Index ndim, const Number* sol, const Number* resid,
                  Index iter, Number norm2_rhs)
  {
    DBG_START_METH("InexactNormalTerminationTester::TestTermination",
                   dbg_verbosity);

    last_iter_ = iter;

    DBG_ASSERT(c_Avc_norm_cauchy_ >= 0.);

    ETerminationTest retval = CONTINUE;

    double norm2_resid = IpBlasDnrm2(ndim, resid, 1);
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "TTNormal: iter = %d ||resid|| = %23.16e ||rhs|| = %23.16e\n",
                   iter, norm2_resid, norm2_rhs);

    if (iter > inexact_normal_max_iter_) {
      IpData().Append_info_string("N+");
      retval = OTHER_SATISFIED;
      return retval;
    }

    //if (Min(norm2_resid/norm2_rhs,norm2_resid) > inexact_normal_tol_) {
    if (norm2_resid/norm2_rhs > inexact_normal_tol_ && norm2_resid > 1e-10) {
      return retval;
    }

    // compute ||c+Av|| for current iterative solution v
    // note that the sol_x and sol_s have the wrong sign
    SmartPtr<const Vector> sol_x;
    SmartPtr<const Vector> sol_s;
    SmartPtr<const Vector> sol_c;
    SmartPtr<const Vector> sol_d;
    GetVectors(ndim, sol, sol_x, sol_s, sol_c, sol_d);

    if (requires_scaling_) {
      SmartPtr<const Vector> scaling_vec = InexCq().curr_scaling_slacks();
      SmartPtr<Vector> tmp = sol_s->MakeNewCopy();
      tmp->ElementWiseMultiply(*scaling_vec);
      sol_s = ConstPtr(tmp);
    }

    SmartPtr<Vector> inf_c = IpCq().curr_c()->MakeNewCopy();
    IpCq().curr_jac_c()->MultVector(-1., *sol_x, 1., *inf_c);
    SmartPtr<const Vector> curr_d_minus_s = IpCq().curr_d_minus_s();
    SmartPtr<Vector> inf_d = curr_d_minus_s->MakeNew();
    inf_d->AddTwoVectors(1., *curr_d_minus_s, 1., *sol_s, 0.);
    IpCq().curr_jac_d()->MultVector(-1., *sol_x, 1., *inf_d);

    Number trial_c_Av_norm = IpCq().CalcNormOfType(NORM_2, *inf_c, *inf_d);

    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "TTNormal: c_Avc_norm_cauchy = %23.16e trial_c_Av_norm = %23.16e\n", c_Avc_norm_cauchy_, trial_c_Av_norm);
    Number BasVal = Max(1., IpCq().curr_primal_infeasibility(NORM_2));

    if (Compare_le(trial_c_Av_norm, c_Avc_norm_cauchy_, BasVal)) {
      retval = OTHER_SATISFIED;
    }

    return retval;
  }

  void
  InexactNormalTerminationTester::Clear()
  {
    DBG_START_METH("InexactNormalTerminationTester::Clear",
                   dbg_verbosity);
  }

} // namespace Ipopt
