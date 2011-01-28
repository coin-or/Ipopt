// Copyright (C) 2009, 2010 International Business Machines and others.
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

#ifdef HAVE_MPI
# include "IpParTripletHelper.hpp"
#endif

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

    // If scaling_factors is NULL, then we are in a parallel setting
    // where the matrix is collected on process 0.  In this case, only
    // process zero collects the scaling_factor entries, but all
    // processes need to call the FillAllValuesFromVector method
    if (scaling_factors) {
      for (Index i=0; i<nx; i++) {
        scaling_factors[i] = 1.;
      }
      scaling_factors += nx;
    }

    SmartPtr<const Vector> scaling_vec = InexCq().curr_scaling_slacks();
#ifdef HAVE_MPI
    ParTripletHelper::FillAllValuesFromVector(ns, *scaling_vec, scaling_factors);
#else
    TripletHelper::FillValuesFromVector(ns, *scaling_vec, scaling_factors);
#endif

    if (scaling_factors) {
      scaling_factors += ns;

      for (Index i=0; i<nc+nd; i++) {
        scaling_factors[i] = 1.;
      }
    }

    return true;
  }

} // namespace Ipopt
