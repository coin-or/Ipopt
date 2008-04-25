// Copyright (C) 2005, 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter                 IBM    2006-01-04
//

#include "IpoptConfig.h"
#include "IpWsmpSolverInterface.hpp"

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
# endif
#endif

/** Prototypes for WSMP's subroutines */
extern "C"
{
  void F77_FUNC(wsetmaxthrds,WSETMAXTHRDS)(const ipfint* NTHREADS);

  void F77_FUNC(wssmp,WSSMP)(const ipfint* N, const ipfint* IA,
                             const ipfint* JA, const double* AVALS,
                             double* DIAG,  ipfint* PERM,
                             ipfint* INVP,  double* B,
                             const ipfint* LDB, const ipfint* NRHS,
                             double* AUX, const ipfint* NAUX,
                             ipfint* MRP, ipfint* IPARM,
                             double* DPARM);
  void F77_FUNC_(wsmp_clear,WSMP_CLEAR)(void);
}

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 3;
#endif

  WsmpSolverInterface::WsmpSolverInterface()
      :
      a_(NULL),
      negevals_(-1),
      initialized_(false),

      PERM_(NULL),
      INVP_(NULL),
      MRP_(NULL)
  {
    DBG_START_METH("WsmpSolverInterface::WsmpSolverInterface()",dbg_verbosity);

    IPARM_ = new ipfint[64];
    DPARM_ = new double[64];
  }

  WsmpSolverInterface::~WsmpSolverInterface()
  {
    DBG_START_METH("WsmpSolverInterface::~WsmpSolverInterface()",
                   dbg_verbosity);

    // Clear WSMP's memory
    F77_FUNC_(wsmp_clear,WSMP_CLEAR)();

    delete[] PERM_;
    delete[] INVP_;
    delete[] MRP_;
    delete[] IPARM_;
    delete[] DPARM_;
    delete[] a_;
  }

  void WsmpSolverInterface::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddLowerBoundedIntegerOption(
      "wsmp_num_threads",
      "Number of threads to be used in WSMP",
      1, 1,
      "This determines on how many processors WSMP is running on.  This option "
      "is only available if Ipopt has been compiled with WSMP.");
    roptions->AddBoundedIntegerOption(
      "wsmp_ordering_option",
      "Determines how ordering is done in WSMP",
      -2, 3, 1,
      "This corresponds to the value of WSSMP's IPARM(16).  This option is "
      "only available if Ipopt has been compiled with WSMP.");
    roptions->AddBoundedNumberOption(
      "wsmp_pivtol",
      "Pivot tolerance for the linear solver WSMP.",
      0.0, true, 1.0, true, 1e-4,
      "A smaller number pivots for sparsity, a larger number pivots for "
      "stability.  This option is only available if Ipopt has been compiled "
      "with WSMP.");
    roptions->AddBoundedNumberOption(
      "wsmp_pivtolmax",
      "Maximum pivot tolerance for the linear solver WSMP.",
      0.0, true, 1.0, true, 1e-1,
      "Ipopt may increase pivtol as high as pivtolmax to get a more accurate "
      "solution to the linear system.  This option is only available if Ipopt "
      "has been compiled with WSMP.");
    roptions->AddBoundedIntegerOption(
      "wsmp_scaling",
      "Determines how the matrix is scaled by WSMP.",
      0, 3, 0,
      "This corresponds to the value of WSSMP's IPARM(10). "
      "This option is only available if Ipopt has been compiled "
      "with WSMP.");
    roptions->AddBoundedNumberOption(
      "wsmp_singularity_threshold",
      "WSMP's singularity threshold.",
      0.0, true, 1.0, true, 1e-18,
      "WSMP's DPARM(10) parameter.  The smaller this value the less likely "
      "a matrix is declared singular.  This option is only available if Ipopt "
      "has been compiled with WSMP.");
    roptions->SetRegisteringCategory("Uncategorized");
    roptions->AddLowerBoundedIntegerOption(
      "wsmp_write_matrix_iteration",
      "Iteration in which the matrices are written to files.",
      -1, -1,
      "If non-negative, this option determines the iteration in which all "
      "matrices given to WSMP are written to files.  This option is only "
      "available if Ipopt has been compiled with WSMP.");
  }

  bool WsmpSolverInterface::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    options.GetIntegerValue("wsmp_num_threads", wsmp_num_threads_, prefix);
    options.GetIntegerValue("wsmp_ordering_option", wsmp_ordering_option_,
                            prefix);
    options.GetNumericValue("wsmp_pivtol", wsmp_pivtol_, prefix);
    if (options.GetNumericValue("wsmp_pivtolmax", wsmp_pivtolmax_, prefix)) {
      ASSERT_EXCEPTION(wsmp_pivtolmax_>=wsmp_pivtol_, OPTION_INVALID,
                       "Option \"wsmp_pivtolmax\": This value must be between "
                       "wsmp_pivtol and 1.");
    }
    else {
      wsmp_pivtolmax_ = Max(wsmp_pivtolmax_, wsmp_pivtol_);
    }
    options.GetNumericValue("wsmp_singularity_threshold",
                            wsmp_singularity_threshold_, prefix);
    options.GetIntegerValue("wsmp_scaling", wsmp_scaling_, prefix);
    options.GetIntegerValue("wsmp_write_matrix_iteration",
                            wsmp_write_matrix_iteration_, prefix);

    // Reset all private data
    dim_=0;
    initialized_=false;
    pivtol_changed_ = false;
    have_symbolic_factorization_ = false;
    delete[] a_;
    a_ = NULL;
    delete[] PERM_;
    PERM_ = NULL;
    delete[] INVP_;
    INVP_ = NULL;
    delete[] MRP_;
    MRP_ = NULL;

    // Set the number of threads
    ipfint NTHREADS = wsmp_num_threads_;
    F77_FUNC(wsetmaxthrds,WSETMAXTHRDS)(&NTHREADS);
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "WSMP will use %d threads.\n", wsmp_num_threads_);

    // Get WSMP's default parameters and set the ones we want differently
    IPARM_[0] = 0;
    IPARM_[1] = 0;
    IPARM_[2] = 0;
    ipfint idmy;
    double ddmy;
    F77_FUNC(wssmp,WSSMP)(&idmy, &idmy, &idmy, &ddmy, &ddmy, &idmy,
                          &idmy, &ddmy, &idmy, &idmy, &ddmy, &idmy,
                          &idmy, IPARM_, DPARM_);
    IPARM_[14] = 0; // no restrictions on pivoting (ignored for
    // Bunch-Kaufman)
    IPARM_[15] = wsmp_ordering_option_; // ordering option
    IPARM_[17] = 1; // use local minimum fill-in ordering
    IPARM_[19] = 1; // for ordering in IP methods?
    IPARM_[30] = 2; // want L D L^T factorization with diagonal
    // pivoting (Bunch/Kaufman)
    //IPARM_[31] = 1; // need D to see where first negative eigenvalue occurs
    //                   if we change this, we need DIAG arguments below!

    IPARM_[10] = 2; // Mark bad pivots

    // Set WSMP's scaling option
    IPARM_[9] = wsmp_scaling_;

    DPARM_[9] = wsmp_singularity_threshold_;

    matrix_file_number_ = 0;

    // Check for SPINLOOPTIME and YIELDLOOPTIME?

    return true;
  }

  ESymSolverStatus WsmpSolverInterface::MultiSolve(
    bool new_matrix,
    const Index* ia,
    const Index* ja,
    Index nrhs,
    double* rhs_vals,
    bool check_NegEVals,
    Index numberOfNegEVals)
  {
    DBG_START_METH("WsmpSolverInterface::MultiSolve",dbg_verbosity);
    DBG_ASSERT(!check_NegEVals || ProvidesInertia());
    DBG_ASSERT(initialized_);

    // check if a factorization has to be done
    if (new_matrix || pivtol_changed_) {
      pivtol_changed_ = false;
      // perform the factorization
      ESymSolverStatus retval;
      retval = Factorization(ia, ja, check_NegEVals, numberOfNegEVals);
      if (retval!=SYMSOLVER_SUCCESS) {
        DBG_PRINT((1, "FACTORIZATION FAILED!\n"));
        return retval;  // Matrix singular or error occurred
      }
    }

    // do the solve
    return Solve(ia, ja, nrhs, rhs_vals);
  }

  double* WsmpSolverInterface::GetValuesArrayPtr()
  {
    DBG_ASSERT(initialized_);
    DBG_ASSERT(a_);
    return a_;
  }

  /** Initialize the local copy of the positions of the nonzero
      elements */
  ESymSolverStatus WsmpSolverInterface::InitializeStructure
  (Index dim, Index nonzeros,
   const Index* ia,
   const Index* ja)
  {
    DBG_START_METH("WsmpSolverInterface::InitializeStructure",dbg_verbosity);
    dim_ = dim;

    // Make space for storing the matrix elements
    delete[] a_;
    a_ = NULL;
    a_ = new double[nonzeros];

    // Do the symbolic facotrization
    ESymSolverStatus retval = SymbolicFactorization(ia, ja);
    if (retval != SYMSOLVER_SUCCESS) {
      return retval;
    }

    initialized_ = true;

    return retval;
  }

  ESymSolverStatus
  WsmpSolverInterface::SymbolicFactorization(
    const Index* ia,
    const Index* ja)
  {
    DBG_START_METH("WsmpSolverInterface::SymbolicFactorization",
                   dbg_verbosity);

    // This is postponed until the first factorization call, since
    // then the values in the matrix are known
    return SYMSOLVER_SUCCESS;
  }

  ESymSolverStatus
  WsmpSolverInterface::InternalSymFact(
    const Index* ia,
    const Index* ja)
  {
    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemSymbolicFactorization().Start();
    }

    // Create space for the permutations
    delete [] PERM_;
    PERM_ = NULL;
    delete [] INVP_;
    INVP_ = NULL;
    delete [] MRP_;
    MRP_ = NULL;
    PERM_ = new ipfint[dim_];
    INVP_ = new ipfint[dim_];
    MRP_ = new ipfint[dim_];

    // Call WSSMP for ordering and symbolic factorization
    ipfint N = dim_;
    ipfint NAUX = 0;
    IPARM_[1] = 1; // ordering
    IPARM_[2] = 2; // symbolic factorization
    ipfint idmy;
    double ddmy;
    F77_FUNC(wssmp,WSSMP)(&N, ia, ja, a_, &ddmy, PERM_, INVP_,
                          &ddmy, &idmy, &idmy, &ddmy, &NAUX, MRP_,
                          IPARM_, DPARM_);

    Index ierror = IPARM_[63];
    if (ierror!=0) {
      if (ierror==-102) {
        Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                       "Error: WSMP is not able to allocate sufficient amount of memory during ordering/symbolic factorization.\n");
      }
      else if (ierror>0) {
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                       "Matrix appears to be singular (with ierror = %d).\n",
                       ierror);
        if (HaveIpData()) {
          IpData().TimingStats().LinearSystemSymbolicFactorization().End();
        }
        return SYMSOLVER_SINGULAR;
      }
      else {
        Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                       "Error in WSMP during ordering/symbolic factorization phase.\n     Error code is %d.\n", ierror);
      }
      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemSymbolicFactorization().End();
      }
      return SYMSOLVER_FATAL_ERROR;
    }
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Predicted memory usage for WSSMP after symbolic factorization IPARM(23)= %d.\n",
                   IPARM_[22]);
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Predicted number of nonzeros in factor for WSSMP after symbolic factorization IPARM(23)= %d.\n",
                   IPARM_[23]);

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemSymbolicFactorization().End();
    }

    return SYMSOLVER_SUCCESS;
  }

  ESymSolverStatus
  WsmpSolverInterface::Factorization(
    const Index* ia,
    const Index* ja,
    bool check_NegEVals,
    Index numberOfNegEVals)
  {
    DBG_START_METH("WsmpSolverInterface::Factorization",dbg_verbosity);

    // If desired, write out the matrix
    Index iter_count = -1;
    if (HaveIpData()) {
      iter_count = IpData().iter_count();
    }
    if (iter_count == wsmp_write_matrix_iteration_) {
      matrix_file_number_++;
      char buf[256];
      sprintf(buf, "wsmp_matrix_%d_%d.dat", iter_count,
              matrix_file_number_);
      Jnlst().Printf(J_SUMMARY, J_LINEAR_ALGEBRA,
                     "Writing WSMP matrix into file %s.\n", buf);
      FILE* fp = fopen(buf, "w");
      fprintf(fp, "%d\n", dim_); // N
      for (Index icol=0; icol<dim_; icol++) {
        fprintf(fp, "%d", ia[icol+1]-ia[icol]); // number of elements for this column
        // Now for each colum we write row indices and values
        for (Index irow=ia[icol]; irow<ia[icol+1]; irow++) {
          fprintf(fp, " %23.16e %d",a_[irow-1],ja[irow-1]);
        }
        fprintf(fp, "\n");
      }
      fclose(fp);
    }

    // Check if we have to do the symbolic factorization and ordering
    // phase yet
    if (!have_symbolic_factorization_) {
      ESymSolverStatus retval = InternalSymFact(ia, ja);
      if (retval != SYMSOLVER_SUCCESS) {
        return retval;
      }
      have_symbolic_factorization_ = true;
    }

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemFactorization().Start();
    }

    // Call WSSMP for numerical factorization
    ipfint N = dim_;
    ipfint NAUX = 0;
    IPARM_[1] = 3; // numerical factorization
    IPARM_[2] = 3; // numerical factorization
    DPARM_[10] = wsmp_pivtol_; // set current pivot tolerance
    ipfint idmy;
    double ddmy;

    F77_FUNC(wssmp,WSSMP)(&N, ia, ja, a_, &ddmy, PERM_, INVP_, &ddmy, &idmy,
                          &idmy, &ddmy, &NAUX, MRP_, IPARM_, DPARM_);

    const Index ierror = IPARM_[63];
    if (ierror > 0) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "WSMP detected that the matrix is singular and encountered %d zero pivots.\n", dim_+1-ierror);
      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemFactorization().End();
      }
      return SYMSOLVER_SINGULAR;
    }
    else if (ierror != 0) {
      if (ierror == -102) {
        Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                       "Error: WSMP is not able to allocate sufficient amount of memory during factorization.\n");
      }
      else {
        Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                       "Error in WSMP during factorization phase.\n     Error code is %d.\n", ierror);
      }
      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemFactorization().End();
      }
      return SYMSOLVER_FATAL_ERROR;
    }
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Memory usage for WSSMP after factorization IPARM(23) = %d\n",
                   IPARM_[22]);
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Number of nonzeros in WSSMP after factorization IPARM(24) = %d\n",
                   IPARM_[23]);

    negevals_ = IPARM_[21]; // Number of negative eigenvalues
    // determined during factorization

    // Check whether the number of negative eigenvalues matches the requested
    // count
    if (check_NegEVals && (numberOfNegEVals!=negevals_)) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Wrong inertia: required are %d, but we got %d.\n",
                     numberOfNegEVals, negevals_);
      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemFactorization().End();
      }
      return SYMSOLVER_WRONG_INERTIA;
    }

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemFactorization().End();
    }
    return SYMSOLVER_SUCCESS;
  }

  ESymSolverStatus WsmpSolverInterface::Solve(
    const Index* ia,
    const Index* ja,
    Index nrhs,
    double *rhs_vals)
  {
    DBG_START_METH("WsmpSolverInterface::Solve",dbg_verbosity);

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemBackSolve().Start();
    }

    // Call WSMP to solve for some right hand sides (including
    // iterative refinement)
    // ToDo: Make iterative refinement an option?
    ipfint N = dim_;
    ipfint LDB = dim_;
    ipfint NRHS = nrhs;
    ipfint NAUX = 0;
    IPARM_[1] = 4; // Forward and Backward Elimintation
    IPARM_[2] = 5; // Iterative refinement
    IPARM_[5] = 1;
    DPARM_[5] = 1e-12;

    ipfint idmy;
    double ddmy;
    F77_FUNC(wssmp,WSSMP)(&N, ia, ja, a_, &ddmy, PERM_, INVP_,
                          rhs_vals, &LDB, &NRHS, &ddmy, &NAUX,
                          MRP_, IPARM_, DPARM_);
    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemBackSolve().End();
    }

    Index ierror = IPARM_[63];
    if (ierror!=0) {
      if (ierror==-102) {
        Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                       "Error: WSMP is not able to allocate sufficient amount of memory during ordering/symbolic factorization.\n");
      }
      else {
        Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                       "Error in WSMP during ordering/symbolic factorization phase.\n     Error code is %d.\n", ierror);
      }
      return SYMSOLVER_FATAL_ERROR;
    }
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Number of iterative refinement steps in WSSMP: %d\n",
                   IPARM_[5]);


    return SYMSOLVER_SUCCESS;
  }

  Index WsmpSolverInterface::NumberOfNegEVals() const
  {
    DBG_START_METH("WsmpSolverInterface::NumberOfNegEVals",dbg_verbosity);
    DBG_ASSERT(negevals_>=0);
    return negevals_;
  }

  bool WsmpSolverInterface::IncreaseQuality()
  {
    DBG_START_METH("WsmpSolverInterface::IncreaseQuality",dbg_verbosity);
    if (wsmp_pivtol_ == wsmp_pivtolmax_) {
      return false;
    }
    pivtol_changed_ = true;

    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Indreasing pivot tolerance for WSMP from %7.2e ",
                   wsmp_pivtol_);
    wsmp_pivtol_ = Min(wsmp_pivtolmax_, pow(wsmp_pivtol_,0.75));
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "to %7.2e.\n",
                   wsmp_pivtol_);
    return true;
  }

  bool WsmpSolverInterface::ProvidesDegeneracyDetection() const
  {
    return true;
  }

  ESymSolverStatus WsmpSolverInterface::
  DetermineDependentRows(const Index* ia, const Index* ja,
                         std::list<Index>& c_deps)
  {
    DBG_START_METH("WsmpSolverInterface::DetermineDependentRows",
                   dbg_verbosity);

    c_deps.clear();

    if (!have_symbolic_factorization_) {
      ESymSolverStatus retval = InternalSymFact(ia, ja);
      if (retval != SYMSOLVER_SUCCESS) {
        return retval;
      }
      have_symbolic_factorization_ = true;
    }

    // Call WSSMP for numerical factorization to detect degenerate
    // rows/columns
    ipfint N = dim_;
    ipfint NAUX = 0;
    IPARM_[1] = 3; // numerical factorization
    IPARM_[2] = 3; // numerical factorization
    DPARM_[10] = wsmp_pivtol_; // set current pivot tolerance
    ipfint idmy;
    double ddmy;

    F77_FUNC(wssmp,WSSMP)(&N, ia, ja, a_, &ddmy, PERM_, INVP_, &ddmy, &idmy,
                          &idmy, &ddmy, &NAUX, MRP_, IPARM_, DPARM_);
    const Index ierror = IPARM_[63];
    if (ierror == 0) {
      Index ndegen = IPARM_[20];
      int ii = 0;
      for (int i=0; i<N; i++) {
        if (MRP_[i] == -1) {
          c_deps.push_back(i);
          ii++;
        }
      }
      DBG_ASSERT(ii == ndegen);
    }
    if (ierror > 0) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "WSMP detected that the matrix is singular and encountered %d zero pivots.\n", dim_+1-ierror);
      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemFactorization().End();
      }
      return SYMSOLVER_SINGULAR;
    }
    else if (ierror != 0) {
      if (ierror == -102) {
        Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                       "Error: WSMP is not able to allocate sufficient amount of memory during factorization.\n");
      }
      else {
        Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                       "Error in WSMP during factorization phase.\n     Error code is %d.\n", ierror);
      }
      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemFactorization().End();
      }
      return SYMSOLVER_FATAL_ERROR;
    }

    return SYMSOLVER_SUCCESS;
  }


} // namespace Ipopt
