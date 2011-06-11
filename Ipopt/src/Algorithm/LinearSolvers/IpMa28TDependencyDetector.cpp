// Copyright (C) 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2007-04-17

#include "IpoptConfig.h"

#ifdef COIN_HAS_HSL
#include "CoinHslConfig.h"
#endif

// if we do not have MA28 in HSL or the linear solver loader, then we want to build the MA28 interface
#if defined(COINHSL_HAS_MA28) || defined(HAVE_LINEARSOLVERLOADER)

#include "IpMa28TDependencyDetector.hpp"

/** Prototypes for MA28's Fortran auxilliary function */
extern "C"
{
  void
  F77_FUNC(ma28part,MA28PART)(ipfint* TASK, ipfint* N, ipfint* M, ipfint* NZ,
                              double* A, ipfint* IROW, ipfint* ICOL,
                              double* PIVTOL, ipfint* FILLFACT, ipfint* IVAR,
                              ipfint* NDEGEN, ipfint* IDEGEN, ipfint* LIW,
                              ipfint* IW, ipfint* LRW, double* RW, ipfint* IERR);
}

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  Ma28TDependencyDetector::Ma28TDependencyDetector()
  {}

  void Ma28TDependencyDetector::
  RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddBoundedNumberOption(
      "ma28_pivtol",
      "Pivot tolerance for linear solver MA28.",
      0.0, true, 1., false, 0.01,
      "This is used when MA28 tries to find the dependent constraints.");
  }

  bool Ma28TDependencyDetector::InitializeImpl(
    const OptionsList& options,
    const std::string& prefix)
  {
    options.GetNumericValue("ma28_pivtol", ma28_pivtol_, prefix);
    return true;
  }

  bool Ma28TDependencyDetector::DetermineDependentRows(
    Index n_rows, Index n_cols, Index n_jac_nz, Number* jac_c_vals,
    Index* jac_c_iRow, Index* jac_c_jCol, std::list<Index>& c_deps)
  {
    DBG_START_METH("Ma28TDependencyDetector::DetermineDependentRows",
                   dbg_verbosity);

    DBG_ASSERT(sizeof(ipfint) == sizeof(Index));

    c_deps.clear();

    // Now comes the interesting part:
    // Call Ma28 to get the dependencies
    ipfint TASK = 0;
    ipfint N = n_cols;
    ipfint M = n_rows;
    ipfint NZ = n_jac_nz;
    double PIVTOL = ma28_pivtol_;
    ipfint FILLFACT = 40;
    ipfint* IVAR;
    ipfint NDEGEN;
    ipfint* IDEGEN;
    ipfint LRW;
    ipfint LIW;
    double ddummy;
    ipfint idummy;
    ipfint IERR;
    // First determine how much work space we need to allocate
    IVAR = new ipfint[N];
    IDEGEN = new ipfint[M];
    F77_FUNC(ma28part,MA28PART)(&TASK, &N, &M, &NZ, &ddummy, jac_c_iRow,
                                jac_c_jCol, &PIVTOL, &FILLFACT, IVAR, &NDEGEN,
                                IDEGEN, &LIW, &idummy, &LRW, &ddummy, &IERR);
    ipfint* IW = new ipfint[LIW];
    double* RW = new double[LRW];

    // Now do the actual factorization and determine dependent constraints
    TASK = 1;
    F77_FUNC(ma28part,MA28PART)(&TASK, &N, &M, &NZ, jac_c_vals, jac_c_iRow,
                                jac_c_jCol, &PIVTOL, &FILLFACT, IVAR, &NDEGEN,
                                IDEGEN, &LIW, IW, &LRW, RW, &IERR);
    delete [] IVAR;
    delete [] IW;
    delete [] RW;
    if (IERR != 0) {
      jnlst_->Printf(J_WARNING, J_INITIALIZATION,
                     "MA28 returns IERR = %d when trying to determine dependent constraints\n", IERR);
      delete [] IDEGEN;
      return false;
    }

    for (Index i=0; i<NDEGEN; i++) {
      c_deps.push_back(IDEGEN[i]-1);
    }

    delete [] IDEGEN;

    return true;
  }

} // namespace Ipopt

#endif /* COINHSL_HAS_MA28 or HAVE_LINEARSOLVERLOADER */
