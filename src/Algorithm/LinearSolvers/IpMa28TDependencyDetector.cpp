// Copyright (C) 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Andreas Waechter            IBM    2007-04-17

#include "IpoptConfig.h"
#include "IpTypes.h"

#ifdef IPOPT_HAS_HSL
#include "CoinHslConfig.h"
#endif

// if we have MA28 in HSL and can compile Fortran code, then we want to build the MA28 interface
#if ((defined(COINHSL_HAS_MA28) && !defined(IPOPT_SINGLE)) || (defined(COINHSL_HAS_MA28S) && defined(IPOPT_SINGLE))) && defined(F77_FUNC)

#include "IpMa28TDependencyDetector.hpp"

/** Prototypes for MA28's Fortran auxiliary function */
extern "C"
{
   void
   F77_FUNC(ma28part, MA28PART)(
      ipindex* TASK,
      ipindex* N,
      ipindex* M,
      ipindex* NZ,
      ipnumber* A,
      ipindex* IROW,
      ipindex* ICOL,
      ipnumber* PIVTOL,
      ipindex* FILLFACT,
      ipindex* IVAR,
      ipindex* NDEGEN,
      ipindex* IDEGEN,
      ipindex* LIW,
      ipindex* IW,
      ipindex* LRW,
      ipnumber* RW,
      ipindex* IERR
   );
}

namespace Ipopt
{
#if IPOPT_VERBOSITY > 0
static const Index dbg_verbosity = 0;
#endif

Ma28TDependencyDetector::Ma28TDependencyDetector()
{ }

void Ma28TDependencyDetector::RegisterOptions(
   SmartPtr<RegisteredOptions> roptions
)
{
   roptions->AddBoundedNumberOption(
      "ma28_pivtol",
      "Pivot tolerance for linear solver MA28.",
      0.0, true,
      1., false,
      0.01);
}

bool Ma28TDependencyDetector::InitializeImpl(
   const OptionsList& options,
   const std::string& prefix
)
{
   options.GetNumericValue("ma28_pivtol", ma28_pivtol_, prefix);
   return true;
}

bool Ma28TDependencyDetector::DetermineDependentRows(
   Index             n_rows,
   Index             n_cols,
   Index             n_jac_nz,
   Number*           jac_c_vals,
   Index*            jac_c_iRow,
   Index*            jac_c_jCol,
   std::list<Index>& c_deps
)
{
   DBG_START_METH("Ma28TDependencyDetector::DetermineDependentRows",
                  dbg_verbosity);

   c_deps.clear();

   // Now comes the interesting part:
   // Call Ma28 to get the dependencies
   Index TASK = 0;
   Index N = n_cols;
   Index M = n_rows;
   Index NZ = n_jac_nz;
   Number PIVTOL = ma28_pivtol_;
   Index FILLFACT = 40;
   Index* IVAR;
   Index NDEGEN;
   Index* IDEGEN;
   Index LRW;
   Index LIW;
   Number ddummy;
   Index idummy;
   Index IERR;
   // First determine how much work space we need to allocate
   IVAR = new Index[N];
   IDEGEN = new Index[M];
   F77_FUNC(ma28part, MA28PART)(&TASK, &N, &M, &NZ, &ddummy, jac_c_iRow, jac_c_jCol, &PIVTOL, &FILLFACT, IVAR, &NDEGEN,
                                IDEGEN, &LIW, &idummy, &LRW, &ddummy, &IERR);
   Index* IW = new Index[LIW];
   Number* RW = new Number[LRW];

   // Now do the actual factorization and determine dependent constraints
   TASK = 1;
   F77_FUNC(ma28part, MA28PART)(&TASK, &N, &M, &NZ, jac_c_vals, jac_c_iRow, jac_c_jCol, &PIVTOL, &FILLFACT, IVAR,
                                &NDEGEN, IDEGEN, &LIW, IW, &LRW, RW, &IERR);
   delete[] IVAR;
   delete[] IW;
   delete[] RW;
   if( IERR != 0 )
   {
      jnlst_->Printf(J_WARNING, J_INITIALIZATION,
                     "MA28 returns IERR = %" IPOPT_INDEX_FORMAT " when trying to determine dependent constraints\n", IERR);
      delete[] IDEGEN;
      return false;
   }

   for( Index i = 0; i < NDEGEN; i++ )
   {
      c_deps.push_back(IDEGEN[i] - 1);
   }

   delete[] IDEGEN;

   return true;
}

} // namespace Ipopt

#endif /* COINHSL_HAS_MA28(s) */
