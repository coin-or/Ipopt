// Copyright (C) 2005, 2010 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
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

#ifdef HAVE_CSTDIO
# include <cstdio>
#else
# ifdef HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
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

#ifdef PARDISO_MATCHING_PREPROCESS
  void smat_reordering_pardiso_wsmp_(const ipfint* N, const ipfint* ia, const ipfint* ja, const double* a_, ipfint* a2, ipfint* ja2,  double* a2_,
                                     ipfint* perm2,  double* scale2, ipfint* tmp2_, ipfint preprocess );
#endif

}


namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 3;
#endif

  WsmpSolverInterface::WsmpSolverInterface()
      :
      a_(NULL),

#ifdef PARDISO_MATCHING_PREPROCESS
      ia2(NULL),
      ja2(NULL),
      a2_(NULL),
      perm2(NULL),
      scale2(NULL),
#endif

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

#ifdef PARDISO_MATCHING_PREPROCESS
    delete[] ia2;
    delete[] ja2;
    delete[] a2_;
    delete[] perm2;
    delete[] scale2;
#endif

    delete[] PERM_;
    delete[] INVP_;
    delete[] MRP_;
    delete[] IPARM_;
    delete[] DPARM_;
    delete[] a_;
  }

  void WsmpSolverInterface::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddIntegerOption(
      "wsmp_num_threads",
      "Number of threads to be used in WSMP",
      1,
      "This determines on how many processors WSMP is running on.  This option "
      "is only available if Ipopt has been compiled with WSMP.");
    roptions->AddBoundedIntegerOption(
      "wsmp_ordering_option",
      "Determines how ordering is done in WSMP (IPARM(16)",
      -2, 3, 1,
      "This corresponds to the value of WSSMP's IPARM(16).  This option is "
      "only available if Ipopt has been compiled with WSMP.");
    roptions->AddBoundedIntegerOption(
      "wsmp_ordering_option2",
      "Determines how ordering is done in WSMP (IPARM(20)",
      0, 3, 1,
      "This corresponds to the value of WSSMP's IPARM(20).  This option is "
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
    roptions->AddStringOption2(
      "wsmp_skip_inertia_check",
      "Always pretent inertia is correct.",
      "no",
      "no", "check inertia",
      "yes", "skip inertia check",
      "Setting this option to \"yes\" essentially disables inertia check. "
      "This option makes the algorithm non-robust and easily fail, but it "
      "might give some insight into the necessity of inertia control.");
    roptions->AddStringOption2(
      "wsmp_no_pivoting",
      "Use the static pivoting option of WSMP.",
      "no",
      "no", "use the regular version",
      "yes", "use static pivoting",
      "Setting this option to \"yes\" means that WSMP instructed not to do "
      "pivoting.  This works only in certain situations (when the Hessian "
      "block is known to be positive definite or when we are using L-BFGS). "
      "It can also lead to a lot of fill-in.");
  }

  bool WsmpSolverInterface::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    options.GetIntegerValue("wsmp_num_threads", wsmp_num_threads_, prefix);
    Index wsmp_ordering_option;
    options.GetIntegerValue("wsmp_ordering_option", wsmp_ordering_option,
                            prefix);
    Index wsmp_ordering_option2;
    options.GetIntegerValue("wsmp_ordering_option2", wsmp_ordering_option2,
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
    options.GetBoolValue("wsmp_skip_inertia_check",
                         skip_inertia_check_, prefix);
    options.GetBoolValue("wsmp_no_pivoting",
                         wsmp_no_pivoting_, prefix);

    // Reset all private data
    dim_=0;
    initialized_=false;
    printed_num_threads_ = false;
    pivtol_changed_ = false;
    have_symbolic_factorization_ = false;
    factorizations_since_recomputed_ordering_ = -1;
    delete[] a_;
    a_ = NULL;
    delete[] PERM_;
    PERM_ = NULL;
    delete[] INVP_;
    INVP_ = NULL;
    delete[] MRP_;
    MRP_ = NULL;

#ifdef PARDISO_MATCHING_PREPROCESS
    delete[] ia2;
    ia2 = NULL;

    delete[] ja2;
    ja2 = NULL;

    delete[] a2_;
    a2_ = NULL;

    delete[] perm2;
    perm2 = NULL;

    delete[] scale2;
    scale2 = NULL;
#endif

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
    IPARM_[15] = wsmp_ordering_option; // ordering option
    IPARM_[17] = 0; // use local minimum fill-in ordering
    IPARM_[19] = wsmp_ordering_option2; // for ordering in IP methods?
    if (wsmp_no_pivoting_) {
      IPARM_[30] = 1; // want L D L^T factorization with diagonal no pivoting
      IPARM_[26] = 1;
    }
    else {
      IPARM_[30] = 2; // want L D L^T factorization with diagonal with pivoting
    }
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

    if (!printed_num_threads_) {
      Jnlst().Printf(J_ITERSUMMARY, J_LINEAR_ALGEBRA,
                     "  -- WSMP is working with %d thread%s.\n", IPARM_[32],
                     IPARM_[32]==1 ? "" : "s");
      printed_num_threads_ = true;
    }
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
    nonzeros_ = nonzeros;

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
    const Index* ja,
    Index numberOfNegEVals)
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

    ipfint N = dim_;

#ifdef PARDISO_MATCHING_PREPROCESS

    delete[] ia2;
    ia2 = NULL;

    delete[] ja2;
    ja2 = NULL;

    delete[] a2_;
    a2_ = NULL;

    delete[] perm2;
    perm2 = NULL;

    delete[] scale2;
    scale2 = NULL;

    ia2    = new ipfint[N+1];
    ja2    = new ipfint[nonzeros_];
    a2_    = new double[nonzeros_];
    perm2  = new ipfint[N];
    scale2 = new double[N];
    ipfint* tmp2_  = new ipfint[N];

    smat_reordering_pardiso_wsmp_(&N, ia, ja, a_, ia2, ja2, a2_, perm2,
                                  scale2, tmp2_, 0);

    delete[] tmp2_;

#endif


    // Call WSSMP for ordering and symbolic factorization
    ipfint NAUX = 0;
    IPARM_[1] = 1; // ordering
    IPARM_[2] = 2; // symbolic factorization
#ifdef PARDISO_MATCHING_PREPROCESS
    IPARM_[9]  =  2; // switch off WSMP's ordering and scaling
    IPARM_[15] = -1; // switch off WSMP's ordering and scaling
    IPARM_[30] =  6; // next step supernode pivoting , since not implemented
    // =2 regular Bunch/Kaufman
    // =1 no pivots
    // =6 limited pivots
    DPARM_[21] = 2e-8; // set pivot perturbation
#endif
    ipfint idmy;
    double ddmy;

    if (wsmp_no_pivoting_) {
      IPARM_[14] = dim_ - numberOfNegEVals; // CHECK
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "Restricting WSMP static pivot sequence with IPARM(15) = %d\n", IPARM_[14]);
    }

    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "Calling WSSMP-1-2 for ordering and symbolic factorization at cpu time %10.3f (wall %10.3f).\n", CpuTime(), WallclockTime());
#ifdef PARDISO_MATCHING_PREPROCESS
    F77_FUNC(wssmp,WSSMP)(&N,  ia2,  ja2,  a2_, &ddmy, PERM_, INVP_,
#else
    F77_FUNC(wssmp,WSSMP)(&N, ia, ja, a_, &ddmy, PERM_, INVP_,
#endif
                          &ddmy, &idmy, &idmy, &ddmy, &NAUX, MRP_,
                          IPARM_, DPARM_);
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "Done with WSSMP-1-2 for ordering and symbolic factorization at cpu time %10.3f (wall %10.3f).\n", CpuTime(), WallclockTime());

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
      Snprintf(buf, 255, "wsmp_matrix_%d_%d.dat", iter_count,
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
      ESymSolverStatus retval = InternalSymFact(ia, ja, numberOfNegEVals);
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

#ifdef PARDISO_MATCHING_PREPROCESS
    {
      ipfint* tmp2_  = new ipfint[N];
      smat_reordering_pardiso_wsmp_ (&N, ia, ja, a_, ia2, ja2, a2_, perm2, scale2, tmp2_, 1);
      delete[] tmp2_;
    }
#endif


    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "Calling WSSMP-3-3 for numerical factorization at cpu time %10.3f (wall %10.3f).\n", CpuTime(), WallclockTime());
#ifdef PARDISO_MATCHING_PREPROCESS
    F77_FUNC(wssmp,WSSMP)(&N, ia2, ja2, a2_, &ddmy, PERM_, INVP_, &ddmy, &idmy,
#else
    F77_FUNC(wssmp,WSSMP)(&N,  ia,  ja,  a_, &ddmy, PERM_, INVP_, &ddmy, &idmy,
#endif
                          &idmy, &ddmy, &NAUX, MRP_, IPARM_, DPARM_);
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "Done with WSSMP-3-3 for numerical factorization at cpu time %10.3f (wall %10.3f).\n", CpuTime(), WallclockTime());

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

    if (factorizations_since_recomputed_ordering_ != -1) {
      factorizations_since_recomputed_ordering_++;
    }

    negevals_ = IPARM_[21]; // Number of negative eigenvalues determined during factorization

    // Check whether the number of negative eigenvalues matches the requested
    // count
    if (check_NegEVals && (numberOfNegEVals!=negevals_)) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Wrong inertia: required are %d, but we got %d.\n",
                     numberOfNegEVals, negevals_);
      if (skip_inertia_check_) {
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                       "  But wsmp_skip_inertia_check is set.  Ignore inertia.\n");
        IpData().Append_info_string("IC ");
        negevals_ = numberOfNegEVals;
      }
      else {
        if (HaveIpData()) {
          IpData().TimingStats().LinearSystemFactorization().End();
        }
        return SYMSOLVER_WRONG_INERTIA;
      }
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

    double ddmy;
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "Calling WSSMP-4-5 for backsolve at cpu time %10.3f (wall %10.3f).\n", CpuTime(), WallclockTime());

#ifdef PARDISO_MATCHING_PREPROCESS
    double* X = new double[nrhs*N];

    // Initialize solution with zero and save right hand side
    for (int i = 0; i < nrhs*N; i++) {
      X[perm2[i]] = scale2[i] * rhs_vals[i];
    }
    F77_FUNC(wssmp,WSSMP)(&N, ia, ja, a_, &ddmy, PERM_, INVP_,
                          X, &LDB, &NRHS, &ddmy, &NAUX,
                          MRP_, IPARM_, DPARM_);
    for (int i = 0; i < N; i++) {
      rhs_vals[i] = scale2[i]*X[perm2[i]];
    }
#else
    F77_FUNC(wssmp,WSSMP)(&N, ia, ja, a_, &ddmy, PERM_, INVP_,
                          rhs_vals, &LDB, &NRHS, &ddmy, &NAUX,
                          MRP_, IPARM_, DPARM_);
#endif

    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "Done with WSSMP-4-5 for backsolve at cpu time %10.3f (wall %10.3f).\n", CpuTime(), WallclockTime());
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

#ifdef PARDISO_MATCHING_PREPROCESS
    delete [] X;
#endif

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

    if (factorizations_since_recomputed_ordering_ == -1 ||
        factorizations_since_recomputed_ordering_ > 2) {
      DPARM_[14] = 1.0;
      pivtol_changed_ = true;
      IpData().Append_info_string("RO ");
      factorizations_since_recomputed_ordering_ = 0;
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Triggering WSMP's recomputation of the ordering for next factorization.\n");
      return true;
    }
    if (wsmp_pivtol_ == wsmp_pivtolmax_) {
      return false;
    }
    pivtol_changed_ = true;

    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Increasing pivot tolerance for WSMP from %7.2e ",
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

    ASSERT_EXCEPTION(!wsmp_no_pivoting_, OPTION_INVALID,
                     "WSMP dependency detection does not work without pivoting.");

    if (!have_symbolic_factorization_) {
      ESymSolverStatus retval = InternalSymFact(ia, ja, 0);
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

#ifdef PARDISO_MATCHING_PREPROCESS
    F77_FUNC(wssmp,WSSMP)(&N, ia2, ja2, a2_, &ddmy, PERM_, INVP_, &ddmy, &idmy,
                          &idmy, &ddmy, &NAUX, MRP_, IPARM_, DPARM_);
#else
    F77_FUNC(wssmp,WSSMP)(&N, ia, ja, a_, &ddmy, PERM_, INVP_, &ddmy, &idmy,
                          &idmy, &ddmy, &NAUX, MRP_, IPARM_, DPARM_);
#endif

    const Index ierror = IPARM_[63];
    if (ierror == 0) {
      int ii = 0;
      for (int i=0; i<N; i++) {
        if (MRP_[i] == -1) {
          c_deps.push_back(i);
          ii++;
        }
      }
      DBG_ASSERT(ii == IPARM_[20]);
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
