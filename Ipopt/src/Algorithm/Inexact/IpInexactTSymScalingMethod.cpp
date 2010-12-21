// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter, Frank E. Curtis         IBM    2009-06-12
//               (based on IpMc19TSymScalingMethod.cpp rev 1204)

#include "IpoptConfig.h"
#include "IpInexactTSymScalingMethod.hpp"
#include "IpTripletHelper.hpp"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif


  bool InexactTSymScalingMethod::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    return true;
  }

  bool InexactTSymScalingMethod::ComputeSymTScalingFactors(Index n,
      Index nnz,
      const ipfint* airn,
      const ipfint* ajcn,
      const double* a,
      double* scaling_factors)
  {
    DBG_START_METH("InexactTSymScalingMethod::ComputeTSymScalingFactors",
                   dbg_verbosity);

    const Index nx = IpData().curr()->x()->Dim();
    const Index ns = IpData().curr()->s()->Dim();
    const Index nc = IpData().curr()->y_c()->Dim();
    const Index nd = IpData().curr()->y_d()->Dim();

    for (Index i=0; i<nx; i++) {
      scaling_factors[i] = 1.;
    }
    scaling_factors += nx;

    SmartPtr<const Vector> scaling_vec = InexCq().curr_scaling_slacks();
    TripletHelper::FillValuesFromVector(ns, *scaling_vec, scaling_factors);
    scaling_factors += ns;

    for (Index i=0; i<nc+nd; i++) {
      scaling_factors[i] = 1.;
    }

    return true;
  }

} // namespace Ipopt
