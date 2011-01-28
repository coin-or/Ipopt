// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter                       IBM    2009-11-13
//               (based on IpInexactTSymScalingMethod.cpp)

#include "IpoptConfig.h"
#include "IpSlackBasedTSymScalingMethod.hpp"
#include "IpTripletHelper.hpp"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif


  bool SlackBasedTSymScalingMethod::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    return true;
  }

  bool SlackBasedTSymScalingMethod::ComputeSymTScalingFactors(Index n,
      Index nnz,
      const ipfint* airn,
      const ipfint* ajcn,
      const double* a,
      double* scaling_factors)
  {
    DBG_START_METH("SlackBasedTSymScalingMethod::ComputeTSymScalingFactors",
                   dbg_verbosity);

    const Index nx = IpData().curr()->x()->Dim();
    const Index ns = IpData().curr()->s()->Dim();
    const Index nc = IpData().curr()->y_c()->Dim();
    const Index nd = IpData().curr()->y_d()->Dim();

    for (Index i=0; i<nx; i++) {
      scaling_factors[i] = 1.;
    }
    scaling_factors += nx;

    SmartPtr<Vector> tmp = IpData().curr()->s()->MakeNew();

    // Lower bounds
    SmartPtr<const Matrix> Pd_L = IpNLP().Pd_L();
    SmartPtr<const Vector> curr_slack_s_L = IpCq().curr_slack_s_L();
    DBG_PRINT_MATRIX(1, "Pd_L", *Pd_L);
    DBG_PRINT_VECTOR(1, "curr_slack_s_L", *curr_slack_s_L);
    Pd_L->MultVector(1., *curr_slack_s_L, 0., *tmp);

    // Upper bounds
    SmartPtr<const Matrix> Pd_U = IpNLP().Pd_U();
    SmartPtr<const Vector> curr_slack_s_U = IpCq().curr_slack_s_U();
    DBG_PRINT_MATRIX(1, "Pd_U", *Pd_U);
    DBG_PRINT_VECTOR(1, "curr_slack_s_U", *curr_slack_s_U);
    Pd_U->MultVector(1., *curr_slack_s_U, 1., *tmp);

    SmartPtr<Vector> tmp2 = tmp->MakeNew();
    const Number slack_scale_max = 1.;
    tmp2->Set(slack_scale_max);
    tmp->ElementWiseMin(*tmp2);

    TripletHelper::FillValuesFromVector(ns, *tmp, scaling_factors);
    scaling_factors += ns;

    for (Index i=0; i<nc+nd; i++) {
      scaling_factors[i] = 1.;
    }

    return true;
  }

} // namespace Ipopt
