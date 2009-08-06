// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter, Sanjeeb Dash          IBM    2009-08-03
//                (based on IpWsmpSolverInterface rev 1535)

#include "IpoptConfig.h"
#include "IpParWsmpSolverInterface.hpp"

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
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
  void F77_FUNC(wsetglobind, WSETGLOBIND)(const ipfint* N,
                                          const ipfint* NUMBERIING,
                                          const ipfint* GLI, ipfint* INFO);
  void F77_FUNC_(pws_xpose_ia,PWS_XPOSE_IA)(const ipfint* N, const ipfint* IA,
      const ipfint* JA, ipfint* NNZ,
      ipfint* IERR);
  void F77_FUNC_(pws_xpose_ja,PWS_XPOSE_JA)(const ipfint* N, const ipfint* IA,
      const ipfint* JA, ipfint* TIA,
      ipfint* TJA, ipfint* IERR);
  void F77_FUNC_(pws_xpose_av,PWS_XPOSE_AV)(const ipfint* N, const ipfint* IA,
      const ipfint* JA,
      const double* AVALS, ipfint* TIA,
      ipfint* TJA, double* TAVALS,
      ipfint* IERR);
  void F77_FUNC_(pwsmp_clear,PWSMP_CLEAR)(void);
}

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 3;
#endif

  ParWsmpSolverInterface::ParWsmpSolverInterface()
      :
      num_local_rows_(-1),
      a_local_(NULL),
      ta_local_(NULL),
      tia_local_(NULL),
      tja_local_(NULL),
      negevals_(-1),
      initialized_(false),

      PERM_(NULL),
      INVP_(NULL),
      MRP_(NULL)
  {
    DBG_START_METH("ParWsmpSolverInterface::ParWsmpSolverInterface()",dbg_verbosity);

    IPARM_ = new ipfint[64];
    DPARM_ = new double[64];
  }

  ParWsmpSolverInterface::~ParWsmpSolverInterface()
  {
    DBG_START_METH("ParWsmpSolverInterface::~ParWsmpSolverInterface()",
                   dbg_verbosity);

    // Clear WSMP's memory
    F77_FUNC_(pwsmp_clear,PWSMP_CLEAR)();

    delete[] PERM_;
    delete[] INVP_;
    delete[] MRP_;
    delete[] IPARM_;
    delete[] DPARM_;
    delete[] a_local_;
    delete[] ta_local_;
    delete[] tia_local_;
    delete[] tja_local_;
  }

  void ParWsmpSolverInterface::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {}

  bool ParWsmpSolverInterface::InitializeImpl(const OptionsList& options,
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

    // Reset all private data
    dim_=0;
    initialized_=false;
    pivtol_changed_ = false;
    have_symbolic_factorization_ = false;
    delete[] a_local_;
    a_local_ = NULL;
    delete[] ta_local_;
    ta_local_ = NULL;
    delete[] tia_local_;
    tia_local_ = NULL;
    delete[] tja_local_;
    tja_local_ = NULL;
    delete[] PERM_;
    PERM_ = NULL;
    delete[] INVP_;
    INVP_ = NULL;
    delete[] MRP_;
    MRP_ = NULL;

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
    IPARM_[17] = 0; // use local minimum fill-in ordering
    IPARM_[19] = wsmp_ordering_option2; // for ordering in IP methods?
    IPARM_[30] = 1; // No pivtoring , since not implemented

    //IPARM_[31] = 1; // need D to see where first negative eigenvalue occurs
    //                   if we change this, we need DIAG arguments below!

    IPARM_[10] = 2; // Mark bad pivots -- WE CAN USE THIS TO DETECT SINGULAR MAIRCES; see also DPARM(21)

    // Set WSMP's scaling option
    IPARM_[9] = wsmp_scaling_;

    DPARM_[9] = wsmp_singularity_threshold_;

    matrix_file_number_ = 0;

    // Check for SPINLOOPTIME and YIELDLOOPTIME?

    return true;
  }

  ESymSolverStatus ParWsmpSolverInterface::MultiSolve(
    bool new_matrix,
    const Index* ia_local,
    const Index* ja_local,
    Index nrhs,
    double* rhs_vals,
    bool check_NegEVals,
    Index numberOfNegEVals)
  {
    DBG_START_METH("ParWsmpSolverInterface::MultiSolve",dbg_verbosity);
    DBG_ASSERT(!check_NegEVals || ProvidesInertia());
    DBG_ASSERT(initialized_);

    // check if a factorization has to be done
    if (new_matrix || pivtol_changed_) {
      pivtol_changed_ = false;
      // perform the factorization
      ESymSolverStatus retval;
      retval = Factorization(ia_local, ja_local, check_NegEVals,
                             numberOfNegEVals);
      if (retval!=SYMSOLVER_SUCCESS) {
        DBG_PRINT((1, "FACTORIZATION FAILED!\n"));
        return retval;  // Matrix singular or error occurred
      }
    }

    // do the solve
    return Solve(ia_local, ja_local, nrhs, rhs_vals);
  }

  double* ParWsmpSolverInterface::GetValuesArrayPtr()
  {
    DBG_ASSERT(initialized_);
    DBG_ASSERT(a_local_);
    return a_local_;
  }

  bool ParWsmpSolverInterface::SetGlobalPos(Index num_rows, const Index* global_pos)
  {
    DBG_ASSERT(num_local_rows_ == -1);
    num_local_rows_ = num_rows;
    const ipfint N = num_rows;
    const ipfint NUMBERING = 1;
    ipfint INFO;
    F77_FUNC(wsetglobind,WSETGLOBIND)(&N, &NUMBERING, global_pos, &INFO);
    if (INFO!=0) {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "wsetglobind returned INFO = %d\n", INFO);
    }
    return (INFO==0);
  }

  /** Initialize the local copy of the positions of the nonzero
      elements */
  ESymSolverStatus ParWsmpSolverInterface::InitializeStructure
  (Index dim, Index nonzeros_local,
   const Index* ia_local,
   const Index* ja_local)
  {
    DBG_START_METH("ParWsmpSolverInterface::InitializeStructure",dbg_verbosity);
    dim_ = dim;

    // Compute the transpose of the system
    DBG_ASSERT(num_local_rows_ != -1);
    const ipfint N = num_local_rows_;
    ipfint NNZ;
    ipfint IERR;
    F77_FUNC_(pws_xpose_ia,PWS_XPOSE_IA)(&N, ia_local, ja_local, &NNZ, &IERR);
    if (IERR!=0) {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "pws_xpose_ia returned IERR = %d\n", IERR);
      return SYMSOLVER_FATAL_ERROR;
    }
    nnz_transpose_local_ = NNZ;
    delete [] tia_local_;
    delete [] tja_local_;
    tia_local_ = NULL;
    tja_local_ = NULL;
    tia_local_ = new Index[nnz_transpose_local_];
    tja_local_ = new Index[nnz_transpose_local_];
    F77_FUNC_(pws_xpose_ja,PWS_XPOSE_JA)(&N, ia_local, ja_local, tia_local_,
                                         tja_local_, &IERR);
    if (IERR!=0) {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "pws_xpose_ja returned IERR = %d\n", IERR);
      return SYMSOLVER_FATAL_ERROR;
    }

    if (Jnlst().ProduceOutput(J_MOREMATRIX, J_LINEAR_ALGEBRA)) {
      const Index* ia_local = tia_local_;
      const Index* ja_local = tja_local_-1;
      Jnlst().StartDistributedOutput();
      int my_rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
      Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA, "On Process %d, WSMP's transpose structure:\n", my_rank);
      for (Index i=0; i<num_local_rows_; i++) {
        Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                       "  tia_local[%5d] = %5d\n", i, ia_local[i]);
        for (Index j=ia_local[i]; j<ia_local[i+1]; j++) {
          Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                         "    tja_local[%5d] = %5d\n", j, ja_local[j]);
        }
      }
      Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                     "  tia_local[%5d] = %5d\n", num_local_rows_, ia_local[num_local_rows_]);
      Jnlst().FinishDistributedOutput();
    }

    // Make space for storing the local matrix elements
    delete[] a_local_;
    a_local_ = NULL;
    a_local_ = new double[nonzeros_local];
    delete[] ta_local_;
    ta_local_ = NULL;
    ta_local_ = new double[nnz_transpose_local_];

    // Do the symbolic facotrization
    ESymSolverStatus retval = SymbolicFactorization(tia_local_,tja_local_);
    if (retval != SYMSOLVER_SUCCESS) {
      return retval;
    }

    initialized_ = true;

    return retval;
  }

  ESymSolverStatus
  ParWsmpSolverInterface::SymbolicFactorization(
    const Index* ia_local,
    const Index* ja_local)
  {
    DBG_START_METH("ParWsmpSolverInterface::SymbolicFactorization",
                   dbg_verbosity);

    // This is postponed until the first factorization call, since
    // then the values in the matrix are known
    return SYMSOLVER_SUCCESS;
  }

  ESymSolverStatus
  ParWsmpSolverInterface::InternalSymFact(
    const Index* ia_local,
    const Index* ja_local,
    Index numberOfNegEVals)
  {
    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemSymbolicFactorization().Start();
    }
    DBG_ASSERT(num_local_rows_ != -1);

    // Create space for the permutations
    delete [] PERM_;
    PERM_ = NULL;
    delete [] INVP_;
    INVP_ = NULL;
    delete [] MRP_;
    MRP_ = NULL;
    PERM_ = new ipfint[dim_];
    INVP_ = new ipfint[dim_];
    MRP_ = new ipfint[num_local_rows_];

//#define TRY
#ifdef TRY

    {
      IPARM_[14] = 0;


      // Call PWSSMP for ordering and symbolic factorization
      ipfint N = num_local_rows_;
      ipfint NAUX = 0;
      IPARM_[1] = 1; // ordering
      IPARM_[2] = 2; // symbolic factorization
      ipfint idmy;
      double ddmy;
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "Calling PWSSMP-1-2 for ordering and symbolic factorization at cpu time %10.3f (wall %10.3f).\n", CpuTime(), WallclockTime());
      F77_FUNC(pwssmp,PWSSMP)(&N, tia_local_, tja_local_, ta_local_, &ddmy, PERM_, INVP_,
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

      Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                     "Predicted number of nonzeros in factor for WSSMP after symbolic factorization IPARM(23)= %d.\n",
                     IPARM_[23]);

    }

#endif

    IPARM_[14] = dim_ - numberOfNegEVals; // CHECK
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "Restricting WSMP static pivot sequence with IPARM(15) = %d\n", IPARM_[14]);

    // Call PWSSMP for ordering and symbolic factorization
    ipfint N = num_local_rows_;
    ipfint NAUX = 0;
    IPARM_[1] = 1; // ordering
    IPARM_[2] = 2; // symbolic factorization
    ipfint idmy;
    double ddmy;
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "Calling PWSSMP-1-2 for ordering and symbolic factorization at cpu time %10.3f (wall %10.3f).\n", CpuTime(), WallclockTime());
    F77_FUNC(pwssmp,PWSSMP)(&N, tia_local_, tja_local_, ta_local_, &ddmy, PERM_, INVP_,
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
  ParWsmpSolverInterface::Factorization(
    const Index* ia_local,
    const Index* ja_local,
    bool check_NegEVals,
    Index numberOfNegEVals)
  {
    DBG_START_METH("ParWsmpSolverInterface::Factorization",dbg_verbosity);

#if 0
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
#endif

    // Transpose the values
    DBG_ASSERT(num_local_rows_ != -1);
    DBG_ASSERT(tia_local_);
    ipfint N = num_local_rows_;
    ipfint IERR;
    F77_FUNC_(pws_xpose_av,PWS_XPOSE_av)(&N, ia_local, ja_local, a_local_,
                                         tia_local_, tja_local_, ta_local_,
                                         &IERR);
    if (IERR!=0) {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "pws_xpose_av returned IERR = %d\n", IERR);
      return SYMSOLVER_FATAL_ERROR;
    }

    // Check if we have to do the symbolic factorization and ordering
    // phase yet
    if (!have_symbolic_factorization_) {
      ESymSolverStatus retval = InternalSymFact(ia_local, ja_local,
                                numberOfNegEVals);
      if (retval != SYMSOLVER_SUCCESS) {
        return retval;
      }
      have_symbolic_factorization_ = true;
    }

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemFactorization().Start();
    }

    // Call WSSMP for numerical factorization
    N = num_local_rows_;
    ipfint NAUX = 0;
    IPARM_[1] = 3; // numerical factorization
    IPARM_[2] = 3; // numerical factorization
    DPARM_[10] = wsmp_pivtol_; // set current pivot tolerance
    ipfint idmy;
    double ddmy;

    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "Calling WSSMP-3-3 for numerical factorization at cpu time %10.3f (wall %10.3f).\n", CpuTime(), WallclockTime());
    F77_FUNC(pwssmp,PWSSMP)(&N, tia_local_, tja_local_, ta_local_,
                            &ddmy, PERM_, INVP_, &ddmy, &idmy,
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

// THE FOLLOWING SHOULD STILL BE CORRECT WITHOUT PIVOTING!
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

  ESymSolverStatus ParWsmpSolverInterface::Solve(
    const Index* ia_local,
    const Index* ja_local,
    Index nrhs,
    double *rhs_vals)
  {
    DBG_START_METH("ParWsmpSolverInterface::Solve",dbg_verbosity);

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemBackSolve().Start();
    }

    // Call WSMP to solve for some right hand sides (including
    // iterative refinement)
    // ToDo: Make iterative refinement an option?
    ipfint N = num_local_rows_;
    ipfint LDB = num_local_rows_;
    ipfint NRHS = nrhs;
    ipfint NAUX = 0;
    IPARM_[1] = 4; // Forward and Backward Elimintation
    IPARM_[2] = 5; // Iterative refinement
    IPARM_[5] = 1;
    DPARM_[5] = 1e-12;

    double ddmy;
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "Calling WSSMP-4-5 for backsolve at cpu time %10.3f (wall %10.3f).\n", CpuTime(), WallclockTime());
    F77_FUNC(pwssmp,PWSSMP)(&N, tia_local_, tja_local_, ta_local_,
                            &ddmy, PERM_, INVP_,
                            rhs_vals, &LDB, &NRHS, &ddmy, &NAUX,
                            MRP_, IPARM_, DPARM_);
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


    return SYMSOLVER_SUCCESS;
  }

  Index ParWsmpSolverInterface::NumberOfNegEVals() const
  {
    DBG_START_METH("ParWsmpSolverInterface::NumberOfNegEVals",dbg_verbosity);
    DBG_ASSERT(false);
    return negevals_;
  }

  bool ParWsmpSolverInterface::IncreaseQuality()
  {
    return false;
#if 0
    DBG_START_METH("ParWsmpSolverInterface::IncreaseQuality",dbg_verbosity);
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
#endif
  }


} // namespace Ipopt
