// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
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
  {}


  bool InexactNormalTerminationTester::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
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
  TestTerminaion(Index ndim, const Number* sol, const Number* resid,
                 Index iter, Number norm2_rhs)
  {
    DBG_START_METH("InexactNormalTerminationTester::TestTerminaion",
                   dbg_verbosity);

    ETerminationTest retval = CONTINUE;

    double norm2_resid = IpBlasDnrm2(ndim, resid, 1);
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "TTNormal: ||resid|| = %23.16e ||rhs|| = %23.16e\n",
                   norm2_resid, norm2_rhs);
    if (Min(norm2_resid/norm2_rhs,norm2_resid) < 1e-6) {
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
