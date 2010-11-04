// Copyright (C) 2010 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter                 IBM    2010-11-01
//               (based on IpWsmpSolverInterface.cpp rev 1729)
//

#include "IpoptConfig.h"
#include "IpParCollectWsmpSolverInterface.hpp"

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

#include "IpMpi.hpp"

/** Prototypes for WSMP's subroutines */
extern "C"
{
  void F77_FUNC(wsetmaxthrds,WSETMAXTHRDS)(const ipfint* NTHREADS);

  void F77_FUNC(pwssmp,PWSSMP)(const ipfint* N, const ipfint* IA,
                               const ipfint* JA, const double* AVALS,
                               double* DIAG,  ipfint* PERM,
                               ipfint* INVP,  double* B,
                               const ipfint* LDB, const ipfint* NRHS,
                               double* AUX, const ipfint* NAUX,
                               ipfint* MRP, ipfint* IPARM,
                               double* DPARM);
  void F77_FUNC_(pwsmp_clear,PWSMP_CLEAR)(void);

#ifdef PARDISO_MATCHING_PREPROCESS
  void smat_reordering_pardiso_wsmp_(const ipfint* N, const ipfint* ia,
                                     const ipfint* ja, const double* a,
                                     ipfint* a2, ipfint* ja2,  double* a2_,
                                     ipfint* perm2,  double* scale2,
                                     ipfint* tmp2_, ipfint preprocess);
#endif



}

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 3;
#endif

  ParCollectWsmpSolverInterface::ParCollectWsmpSolverInterface()
      :
      a_(NULL),
      negevals_(-1),
      initialized_(false),

#ifdef PARDISO_MATCHING_PREPROCESS
      ia2_(NULL),
      ja2_(NULL),
      a2_(NULL),
      perm2_(NULL),
      scale2_(NULL),
#endif


      PERM_(NULL),
      INVP_(NULL),
      MRP_(NULL)

  {
    DBG_START_METH("ParCollectWsmpSolverInterface::ParCollectWsmpSolverInterface()",dbg_verbosity);

    IPARM_ = new ipfint[64];
    DPARM_ = new double[64];
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_);
  }

  ParCollectWsmpSolverInterface::~ParCollectWsmpSolverInterface()
  {
    DBG_START_METH("ParCollectWsmpSolverInterface::~ParCollectWsmpSolverInterface()",
                   dbg_verbosity);

    // Clear WSMP's memory
    F77_FUNC_(pwsmp_clear,PWSMP_CLEAR)();

#ifdef PARDISO_MATCHING_PREPROCESS
    delete[] ia2_;
    delete[] ja2_;
    delete[] a2_;
    delete[] perm2_;
    delete[] scale2_;
#endif


    delete[] PERM_;
    delete[] INVP_;
    delete[] MRP_;
    delete[] IPARM_;
    delete[] DPARM_;
    delete[] a_;
  }

  void ParCollectWsmpSolverInterface::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {}

  bool ParCollectWsmpSolverInterface::InitializeImpl(const OptionsList& options,
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
    delete[] ia2_;
    ia2_ = NULL;

    delete[] ja2_;
    ja2_ = NULL;

    delete[] a2_;
    a2_ = NULL;

    delete[] perm2_;
    perm2_ = NULL;

    delete[] scale2_;
    scale2_ = NULL;
#endif


    // Set the number of threads
    if (wsmp_num_threads_==0) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "WSMP uses its defaults number of threads.\n");
    }
    else {
      ipfint NTHREADS = wsmp_num_threads_;
      F77_FUNC(wsetmaxthrds,WSETMAXTHRDS)(&NTHREADS);
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "WSMP will use %d threads.\n", wsmp_num_threads_);
    }

    // Get WSMP's default parameters and set the ones we want differently
    IPARM_[0] = 0;
    IPARM_[1] = 0;
    IPARM_[2] = 0;
    ipfint idmy;
    double ddmy;
    F77_FUNC(pwssmp,PWSSMP)(&idmy, &idmy, &idmy, &ddmy, &ddmy, &idmy,
                            &idmy, &ddmy, &idmy, &idmy, &ddmy, &idmy,
                            &idmy, IPARM_, DPARM_);
    IPARM_[15] = wsmp_ordering_option; // ordering option

#ifdef PARDISO_MATCHING_PREPROCESS
    //IPARM_[9]  =  2; // next step switch off WSMP's ordering and scaling
    //IPARM_[15] = -1; // next step switch off WSMP's ordering and scaling
    //IPARM_[30] =  6; // next step supernode pivoting , since not implemented
#endif

    IPARM_[17] = 0; // use local minimum fill-in ordering
    IPARM_[19] = wsmp_ordering_option2; // for ordering in IP methods?
    IPARM_[30] = 1; // No pivoting , since not implemented
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

  ESymSolverStatus ParCollectWsmpSolverInterface::MultiSolve(
    bool new_matrix,
    const Index* ia,
    const Index* ja,
    Index nrhs,
    double* rhs_vals,
    bool check_NegEVals,
    Index numberOfNegEVals)
  {
    DBG_START_METH("ParCollectWsmpSolverInterface::MultiSolve",dbg_verbosity);
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

  double* ParCollectWsmpSolverInterface::GetValuesArrayPtr()
  {
    DBG_ASSERT(initialized_);
    return a_;
  }

  /** Initialize the local copy of the positions of the nonzero
      elements */
  ESymSolverStatus ParCollectWsmpSolverInterface::InitializeStructure
  (Index dim, Index nonzeros,
   const Index* ia,
   const Index* ja)
  {
    DBG_START_METH("ParCollectWsmpSolverInterface::InitializeStructure",dbg_verbosity);
    dim_ = dim;
    nonzeros_ = nonzeros;

    // Make space for storing the matrix elements
    delete[] a_;
    a_ = NULL;
    if (my_rank_==0) {
      a_ = new double[nonzeros];
    }

    // Do the symbolic facotrization
    ESymSolverStatus retval = SymbolicFactorization(ia, ja);
    if (retval != SYMSOLVER_SUCCESS) {
      return retval;
    }

    initialized_ = true;

    return retval;
  }

  ESymSolverStatus
  ParCollectWsmpSolverInterface::SymbolicFactorization(
    const Index* ia,
    const Index* ja)
  {
    DBG_START_METH("ParCollectWsmpSolverInterface::SymbolicFactorization",
                   dbg_verbosity);

    // This is postponed until the first factorization call, since
    // then the values in the matrix are known
    return SYMSOLVER_SUCCESS;
  }

  ESymSolverStatus
  ParCollectWsmpSolverInterface::InternalSymFact(
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
    if (my_rank_==0) {
      MRP_ = new ipfint[dim_];
    }


    // Call WSSMP for ordering and symbolic factorization
    ipfint N = 0;
    const ipfint one = 1;
    const ipfint* ia_wsmp = &one;
    if (my_rank_==0) {
      N = dim_;
      ia_wsmp = ia;

#ifdef PARDISO_MATCHING_PREPROCESS
      delete[] ia2_;
      ia2_ = NULL;

      delete[] ja2_;
      ja2_ = NULL;

      delete[] a2_;
      a2_ = NULL;

      delete[] perm2_;
      perm2_ = NULL;

      delete[] scale2_;
      scale2_ = NULL;

      ia2_    = new ipfint[N+1];
      ja2_    = new ipfint[nonzeros_];
      a2_     = new double[nonzeros_];
      perm2_  = new ipfint[N];
      scale2_ = new double[N];
      ipfint* tmp2_  = new ipfint[N];

      smat_reordering_pardiso_wsmp_(&N, ia, ja, a_, ia2_, ja2_, a2_,
                                    perm2_, scale2_, tmp2_, 0);

      delete[] tmp2_;
      ia_wsmp = ia2_;

#endif

    }
    ipfint NAUX = 0;
    IPARM_[1] = 1; // ordering
    IPARM_[2] = 2; // symbolic factorization
    ipfint idmy;
    double ddmy;

    IPARM_[14] = dim_ - numberOfNegEVals; // CHECK
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "Restricting WSMP static pivot sequence with IPARM(15) = %d\n", IPARM_[14]);

    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "Calling PWSSMP-1-2 for ordering and symbolic factorization at cpu time %10.3f (wall %10.3f).\n", CpuTime(), WallclockTime());

#ifdef PARDISO_MATCHING_PREPROCESS
    F77_FUNC(pwssmp,PWSSMP)(&N, ia_wsmp, ja2_, a2_, &ddmy, PERM_, INVP_,
#else
    F77_FUNC(pwssmp,PWSSMP)(&N, ia_wsmp, ja, a_, &ddmy, PERM_, INVP_,
#endif
                            &ddmy, &idmy, &idmy, &ddmy, &NAUX, MRP_,
                            IPARM_, DPARM_);

    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "Done with PWSSMP-1-2 for ordering and symbolic factorization at cpu time %10.3f (wall %10.3f).\n", CpuTime(), WallclockTime());

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
                   "Predicted memory usage for PWSSMP after symbolic factorization IPARM(23)= %d.\n",
                   IPARM_[22]);
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Predicted number of nonzeros in factor for PWSSMP after symbolic factorization IPARM(23)= %d.\n",
                   IPARM_[23]);

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemSymbolicFactorization().End();
    }

    return SYMSOLVER_SUCCESS;
  }

  ESymSolverStatus
  ParCollectWsmpSolverInterface::Factorization(
    const Index* ia,
    const Index* ja,
    bool check_NegEVals,
    Index numberOfNegEVals)
  {
    DBG_START_METH("ParCollectWsmpSolverInterface::Factorization",dbg_verbosity);

    // If desired, write out the matrix
    Index iter_count = -1;
    if (HaveIpData()) {
      iter_count = IpData().iter_count();
    }
    if (iter_count == wsmp_write_matrix_iteration_) {
      int mpi_size;
      MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
      ASSERT_EXCEPTION(mpi_size == 1, OPTION_INVALID,
                       "wsmp_write_matrix_iteration is set, but we have more than 1 MPI process. Ouptut of WSMP matrix only implemented for 1 process.\n");
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
    ipfint N = 0;
    const ipfint one = 1;
    const ipfint* ia_wsmp = &one;
    if (my_rank_==0) {
      N = dim_;
#ifdef PARDISO_MATCHING_PREPROCESS
      ia_wsmp = ia2_;
#else
      ia_wsmp = ia;
#endif
    }
    ipfint NAUX = 0;
    IPARM_[1] = 3; // numerical factorization
    IPARM_[2] = 3; // numerical factorization
    DPARM_[10] = wsmp_pivtol_; // set current pivot tolerance
    ipfint idmy;
    double ddmy;

    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "Calling PWSSMP-3-3 for numerical factorization at cpu time %10.3f (wall %10.3f).\n", CpuTime(), WallclockTime());
#ifdef PARDISO_MATCHING_PREPROCESS
    F77_FUNC(pwssmp,PWSSMP)(&N, ia_wsmp, ja2_, a2_, &ddmy, PERM_, INVP_, &ddmy, &idmy,
#else
    F77_FUNC(pwssmp,PWSSMP)(&N, ia_wsmp, ja, a_, &ddmy, PERM_, INVP_, &ddmy, &idmy,
#endif
                            &idmy, &ddmy, &NAUX, MRP_, IPARM_, DPARM_);

    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "Done with PWSSMP-3-3 for numerical factorization at cpu time %10.3f (wall %10.3f).\n", CpuTime(), WallclockTime());

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
                   "Memory usage for PWSSMP after factorization IPARM(23) = %d\n",
                   IPARM_[22]);
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Number of nonzeros in PWSSMP after factorization IPARM(24) = %d\n",
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

  ESymSolverStatus ParCollectWsmpSolverInterface::Solve(
    const Index* ia,
    const Index* ja,
    Index nrhs,
    double *rhs_vals)
  {
    DBG_START_METH("ParCollectWsmpSolverInterface::Solve",dbg_verbosity);

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemBackSolve().Start();
    }

    // Call WSMP to solve for some right hand sides (including
    // iterative refinement)
    // ToDo: Make iterative refinement an option?
    ipfint N = 0;
    ipfint LDB = 0;
    const ipfint one = 1;
    const ipfint* ia_wsmp = &one;
    if (my_rank_==0) {
      N = dim_;
      LDB = dim_;
#ifdef PARDISO_MATCHING_PREPROCESS
      ia_wsmp = ia2_;
#else
      ia_wsmp = ia;
#endif
    }
    ipfint NRHS = nrhs;
    ipfint NAUX = 0;
    IPARM_[1] = 4; // Forward and Backward Elimintation
    IPARM_[2] = 5; // Iterative refinement
    IPARM_[5] = 1;
    DPARM_[5] = 1e-12;

    double ddmy;
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "Calling PWSSMP-4-5 for backsolve at cpu time %10.3f (wall %10.3f).\n", CpuTime(), WallclockTime());

#ifdef PARDISO_MATCHING_PREPROCESS
    double* X;
    if (my_rank_==0) {
      X = new double[nrhs*N];

      // Initialize solution with zero and save right hand side
      for (int i = 0; i < nrhs*N; i++) {
        X[perm2_[i]] = scale2_[i] * rhs_vals[i];
      }
    }
    F77_FUNC(pwssmp,PWSSMP)(&N, ia_wsmp, ja2_, a2_, &ddmy, PERM_, INVP_,
                            X, &LDB, &NRHS, &ddmy, &NAUX,
                            MRP_, IPARM_, DPARM_);
    if (my_rank_==0) {
      for (int i = 0; i < N; i++) {
        rhs_vals[i] = scale2_[i]*X[perm2_[i]];
      }
    }
#else
    F77_FUNC(pwssmp,PWSSMP)(&N, ia_wsmp, ja, a_, &ddmy, PERM_, INVP_,
                            rhs_vals, &LDB, &NRHS, &ddmy, &NAUX,
                            MRP_, IPARM_, DPARM_);
#endif

    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "Done with PWSSMP-4-5 for backsolve at cpu time %10.3f (wall %10.3f).\n", CpuTime(), WallclockTime());
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
                   "Number of iterative refinement steps in PWSSMP: %d\n",
                   IPARM_[5]);

#ifdef PARDISO_MATCHING_PREPROCESS
    delete [] X;
#endif

    return SYMSOLVER_SUCCESS;
  }

  Index ParCollectWsmpSolverInterface::NumberOfNegEVals() const
  {
    DBG_START_METH("ParCollectWsmpSolverInterface::NumberOfNegEVals",dbg_verbosity);
    DBG_ASSERT(negevals_>=0);
    return negevals_;
  }

  bool ParCollectWsmpSolverInterface::IncreaseQuality()
  {
// TODO: BEFORE INCREASING TOLERANCE, TRY TO REDO THE ORDERING. - USE DPARM(15) !!!
//       ALSO: DECREASE PIVTOL AGAIN LATER?
    DBG_START_METH("ParCollectWsmpSolverInterface::IncreaseQuality",dbg_verbosity);

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

} // namespace Ipopt
