// Copyright (C) 2008, 2011 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter             IBM    2008-09-19
//
//           based on IpPardisoSolverInterface.cpp rev 1219

#include "IpoptConfig.h"
#include "IpIterativePardisoSolverInterface.hpp"
#include "IpBlas.hpp"
# include <math.h>

#ifdef HAVE_CSTDIO
# include <cstdio>
#else
# ifdef HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif

#ifdef HAVE_CSTDLIB
# include <cstdlib>
#else
# ifdef HAVE_STDLIB_H
#  include <stdlib.h>
# else
#  error "don't have header file for stdlib"
# endif
#endif

#ifdef HAVE_CSTRING
# include <cstring>
#else
# ifdef HAVE_STRING_H
#  include <string.h>
# else
#  error "don't have header file for string"
# endif
#endif


Ipopt::IterativeSolverTerminationTester* global_tester_ptr_;
Ipopt::IterativeSolverTerminationTester::ETerminationTest test_result_;
extern "C"
{
  int IpoptTerminationTest(int n, double* sol, double* resid, int iter, double norm2_rhs) {
    fflush(stdout);
    fflush(stderr);
    // only do the termination test for PD system
    test_result_ = global_tester_ptr_->TestTermination(n, sol, resid, iter, norm2_rhs);
    global_tester_ptr_->GetJnlst().Printf(Ipopt::J_DETAILED, Ipopt::J_LINEAR_ALGEBRA,
                                          "Termination Tester Result = %d.\n",
                                          test_result_);
    switch (test_result_) {
    case Ipopt::IterativeSolverTerminationTester::CONTINUE:
      return false;
      break;
    default:
      return true;
      break;
    }
  }

  // The following global function pointer is defined in the Pardiso library
  void SetIpoptCallbackFunction(int (*IpoptFunction)(int n, double* x,  double* r, int k, double b));
}

/** Prototypes for Pardiso's subroutines */
extern "C"
{
  void F77_FUNC(pardisoinit,PARDISOINIT)(void* PT, const ipfint* MTYPE,
                                         const ipfint* SOLVER,
                                         ipfint* IPARM,
                                         double* DPARM,
                                         ipfint* ERROR);
  void F77_FUNC(pardiso,PARDISO)(void** PT, const ipfint* MAXFCT,
                                 const ipfint* MNUM, const ipfint* MTYPE,
                                 const ipfint* PHASE, const ipfint* N,
                                 const double* A, const ipfint* IA,
                                 const ipfint* JA, const ipfint* PERM,
                                 const ipfint* NRHS, ipfint* IPARM,
                                 const ipfint* MSGLVL, double* B, double* X,
                                 ipfint* ERROR, double* DPARM);
}

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  IterativePardisoSolverInterface::
  IterativePardisoSolverInterface(IterativeSolverTerminationTester& normal_tester,
                                  IterativeSolverTerminationTester& pd_tester)
      :
      a_(NULL),
      negevals_(-1),
      initialized_(false),

      MAXFCT_(1),
      MNUM_(1),
      MTYPE_(-2),
      MSGLVL_(0),
      debug_last_iter_(-1),
      normal_tester_(&normal_tester),
      pd_tester_(&pd_tester)
  {
    DBG_START_METH("IterativePardisoSolverInterface::IterativePardisoSolverInterface()",dbg_verbosity);

    PT_ = new void*[64];
    IPARM_ = new ipfint[64];
    DPARM_ = new double[64];
  }

  IterativePardisoSolverInterface::~IterativePardisoSolverInterface()
  {
    DBG_START_METH("IterativePardisoSolverInterface::~IterativePardisoSolverInterface()",
                   dbg_verbosity);

    // Tell Pardiso to release all memory
    if (initialized_) {
      ipfint PHASE = -1;
      ipfint N = dim_;
      ipfint NRHS = 0;
      ipfint ERROR;
      ipfint idmy;
      double ddmy;
      F77_FUNC(pardiso,PARDISO)(PT_, &MAXFCT_, &MNUM_, &MTYPE_, &PHASE, &N,
                                &ddmy, &idmy, &idmy, &idmy, &NRHS, IPARM_,
                                &MSGLVL_, &ddmy, &ddmy, &ERROR, DPARM_);
      DBG_ASSERT(ERROR==0);
    }

    delete[] PT_;
    delete[] IPARM_;
    delete[] DPARM_;
    delete[] a_;
  }

  void IterativePardisoSolverInterface::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {}

  bool IterativePardisoSolverInterface::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
#ifdef HAVE_PARDISO_OLDINTERFACE
    THROW_EXCEPTION(OPTION_INVALID, "The inexact version works only with a new version of Pardiso (at least 4.0)");
#endif

    Index enum_int;
    options.GetEnumValue("pardiso_matching_strategy", enum_int, prefix);
    match_strat_ = PardisoMatchingStrategy(enum_int);
    options.GetBoolValue("pardiso_redo_symbolic_fact_only_if_inertia_wrong",
                         pardiso_redo_symbolic_fact_only_if_inertia_wrong_,
                         prefix);
    options.GetBoolValue("pardiso_repeated_perturbation_means_singular",
                         pardiso_repeated_perturbation_means_singular_,
                         prefix);
    //Index pardiso_out_of_core_power;
    //options.GetIntegerValue("pardiso_out_of_core_power",
    //                        pardiso_out_of_core_power, prefix);
    options.GetBoolValue("pardiso_skip_inertia_check",
                         skip_inertia_check_, prefix);
    int max_iterref_steps;
    options.GetIntegerValue("pardiso_max_iterative_refinement_steps", max_iterref_steps, prefix);

    // PD system
    options.GetIntegerValue("pardiso_max_iter", pardiso_max_iter_, prefix);
    options.GetNumericValue("pardiso_iter_relative_tol",
                            pardiso_iter_relative_tol_, prefix);
    options.GetIntegerValue("pardiso_iter_coarse_size",
                            pardiso_iter_coarse_size_, prefix);
    options.GetIntegerValue("pardiso_iter_max_levels",
                            pardiso_iter_max_levels_, prefix);
    options.GetNumericValue("pardiso_iter_dropping_factor",
                            pardiso_iter_dropping_factor_, prefix);
    options.GetNumericValue("pardiso_iter_dropping_schur",
                            pardiso_iter_dropping_schur_, prefix);
    options.GetIntegerValue("pardiso_iter_max_row_fill",
                            pardiso_iter_max_row_fill_, prefix);
    options.GetNumericValue("pardiso_iter_inverse_norm_factor",
                            pardiso_iter_inverse_norm_factor_, prefix);
    // Normal system
    options.GetIntegerValue("pardiso_max_iter", normal_pardiso_max_iter_,
                            prefix+"normal.");
    options.GetNumericValue("pardiso_iter_relative_tol",
                            normal_pardiso_iter_relative_tol_,
                            prefix+"normal.");
    options.GetIntegerValue("pardiso_iter_coarse_size",
                            normal_pardiso_iter_coarse_size_,
                            prefix+"normal.");
    options.GetIntegerValue("pardiso_iter_max_levels",
                            normal_pardiso_iter_max_levels_,
                            prefix+"normal.");
    options.GetNumericValue("pardiso_iter_dropping_factor",
                            normal_pardiso_iter_dropping_factor_,
                            prefix+"normal.");
    options.GetNumericValue("pardiso_iter_dropping_schur",
                            normal_pardiso_iter_dropping_schur_,
                            prefix+"normal.");
    options.GetIntegerValue("pardiso_iter_max_row_fill",
                            normal_pardiso_iter_max_row_fill_,
                            prefix+"normal.");
    options.GetNumericValue("pardiso_iter_inverse_norm_factor",
                            normal_pardiso_iter_inverse_norm_factor_,
                            prefix+"normal.");

    int pardiso_msglvl;
    options.GetIntegerValue("pardiso_msglvl", pardiso_msglvl, prefix);
    int order;
    options.GetEnumValue("pardiso_order", order, prefix);
    options.GetIntegerValue("pardiso_max_droptol_corrections",
                            pardiso_max_droptol_corrections_, prefix);

    // Number value = 0.0;

    // Tell Pardiso to release all memory if it had been used before
    if (initialized_) {
      ipfint PHASE = -1;
      ipfint N = dim_;
      ipfint NRHS = 0;
      ipfint ERROR;
      ipfint idmy;
      double ddmy;
      F77_FUNC(pardiso,PARDISO)(PT_, &MAXFCT_, &MNUM_, &MTYPE_, &PHASE, &N,
                                &ddmy, &idmy, &idmy, &idmy, &NRHS, IPARM_,
                                &MSGLVL_, &ddmy, &ddmy, &ERROR, DPARM_) ;
      DBG_ASSERT(ERROR==0);
    }

    // Reset all private data
    dim_=0;
    nonzeros_=0;
    have_symbolic_factorization_=false;
    initialized_=false;
    delete[] a_;
    a_ = NULL;

    // Call Pardiso's initialization routine
    IPARM_[0] = 0;  // Tell it to fill IPARM with default values(?)
    ipfint ERROR = 0;
    ipfint SOLVER = 1; // initialze only direct solver
    F77_FUNC(pardisoinit,PARDISOINIT)(PT_, &MTYPE_, &SOLVER,
                                      IPARM_, DPARM_, &ERROR);

    // Set some parameters for Pardiso
    IPARM_[0] = 1;  // Don't use the default values

#if defined(HAVE_PARDISO_PARALLEL) || ! defined(HAVE_PARDISO)
    // Obtain the numbers of processors from the value of OMP_NUM_THREADS
    char    *var = getenv("OMP_NUM_THREADS");
    int      num_procs = 1;
    if (var != NULL) {
      sscanf( var, "%d", &num_procs );
      if (num_procs < 1) {
        Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                       "Invalid value for OMP_NUM_THREADS (\"%s\").\n", var);
        return false;
      }
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Using environment OMP_NUM_THREADS = %d as the number of processors.\n", num_procs);
    }
#ifdef HAVE_PARDISO
    // If we run Pardiso through the linear solver loader,
    // we do not know whether it is the parallel version, so we do not report an error if OMP_NUM_THREADS is not set.
    else {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "You need to set environment variable OMP_NUM_THREADS to the number of processors used in Pardiso (e.g., 1).\n\n");
      return false;
    }
#endif
    IPARM_[2] = num_procs;  // Set the number of processors
#else

    IPARM_[2] = 1;
#endif

    IPARM_[1] = order;
    IPARM_[5] = 1;  // Overwrite right-hand side
    IPARM_[7] = max_iterref_steps;

    // Options suggested by Olaf Schenk
    IPARM_[9] = 12;
    IPARM_[10] = 2; // Results in better scaling
    // Matching information:  IPARM_[12] = 1 seems ok, but results in a
    // large number of pivot perturbation
    // Matching information:  IPARM_[12] = 2 robust,  but more  expensive method
    IPARM_[12] = (int)match_strat_;
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Pardiso matching strategy (IPARM(13)): %d\n", IPARM_[12]);

    IPARM_[20] = 3; // Results in better accuracy
    IPARM_[23] = 1; // parallel fac
    IPARM_[24] = 1; // parallel solve
    IPARM_[28] = 0; // 32-bit factorization
    IPARM_[29] = 1; // we need this for IPOPT interface

    IPARM_[31] = 1 ;  // iterative solver
    MSGLVL_ = pardiso_msglvl;

    pardiso_iter_dropping_factor_used_ = pardiso_iter_dropping_factor_;
    pardiso_iter_dropping_schur_used_ = pardiso_iter_dropping_schur_;
    normal_pardiso_iter_dropping_factor_used_ = normal_pardiso_iter_dropping_factor_;
    normal_pardiso_iter_dropping_schur_used_ = normal_pardiso_iter_dropping_schur_;

    // TODO Make option
    decr_factor_ = 1./3.;

    // Option for the out of core variant
    // IPARM_[49] = pardiso_out_of_core_power;

    SetIpoptCallbackFunction(&IpoptTerminationTest);

    bool retval = normal_tester_->Initialize(Jnlst(), IpNLP(), IpData(),
                  IpCq(), options, prefix);
    if (retval) {
      retval = pd_tester_->Initialize(Jnlst(), IpNLP(), IpData(),
                                      IpCq(), options, prefix);
    }

    return retval;
  }

  ESymSolverStatus IterativePardisoSolverInterface::MultiSolve(bool new_matrix,
      const Index* ia,
      const Index* ja,
      Index nrhs,
      double* rhs_vals,
      bool check_NegEVals,
      Index numberOfNegEVals)
  {
    DBG_START_METH("IterativePardisoSolverInterface::MultiSolve",dbg_verbosity);
    DBG_ASSERT(!check_NegEVals || ProvidesInertia());
    DBG_ASSERT(initialized_);

    // check if a factorization has to be done
    if (new_matrix) {
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

  double* IterativePardisoSolverInterface::GetValuesArrayPtr()
  {
    DBG_ASSERT(initialized_);
    DBG_ASSERT(a_);
    return a_;
  }

  /** Initialize the local copy of the positions of the nonzero
      elements */
  ESymSolverStatus IterativePardisoSolverInterface::InitializeStructure
  (Index dim, Index nonzeros,
   const Index* ia,
   const Index* ja)
  {
    DBG_START_METH("IterativePardisoSolverInterface::InitializeStructure",dbg_verbosity);
    dim_ = dim;
    nonzeros_ = nonzeros;

    // Make space for storing the matrix elements
    delete[] a_;
    a_ = NULL;
    a_ = new double[nonzeros_];

    // Do the symbolic facotrization
    ESymSolverStatus retval = SymbolicFactorization(ia, ja);
    if (retval != SYMSOLVER_SUCCESS) {
      return retval;
    }

    initialized_ = true;

    return retval;
  }

  ESymSolverStatus
  IterativePardisoSolverInterface::SymbolicFactorization(const Index* ia,
      const Index* ja)
  {
    DBG_START_METH("IterativePardisoSolverInterface::SymbolicFactorization",
                   dbg_verbosity);

    // Since Pardiso requires the values of the nonzeros of the matrix
    // for an efficient symbolic factorization, we postpone that task
    // until the first call of Factorize.  All we do here is to reset
    // the flag (in case this interface is called for a matrix with a
    // new structure).

    have_symbolic_factorization_ = false;

    return SYMSOLVER_SUCCESS;
  }

  static void
  write_iajaa_matrix (int     N,
                      const Index*  ia,
                      const Index*  ja,
                      double*      a_,
                      double*      rhs_vals,
                      int        iter_cnt,
                      int        sol_cnt)
  {
    if (getenv ("IPOPT_WRITE_MAT")) {
      /* Write header */
      FILE    *mat_file;
      char     mat_name[128];
      char     mat_pref[32];

      ipfint   NNZ = ia[N]-1;
      ipfint   i;

      if (getenv ("IPOPT_WRITE_PREFIX"))
        strcpy (mat_pref, getenv ("IPOPT_WRITE_PREFIX"));
      else
        strcpy (mat_pref, "mat-ipopt");

      Snprintf (mat_name, 127, "%s_%03d-%02d.iajaa",
                mat_pref, iter_cnt, sol_cnt);

      // Open and write matrix file.
      mat_file = fopen (mat_name, "w");

      fprintf (mat_file, "%d\n", N);
      fprintf (mat_file, "%d\n", NNZ);

      for (i = 0; i < N+1; i++)
        fprintf (mat_file, "%d\n", ia[i]);
      for (i = 0; i < NNZ; i++)
        fprintf (mat_file, "%d\n", ja[i]);
      for (i = 0; i < NNZ; i++)
        fprintf (mat_file, "%32.24e\n", a_[i]);

      /* Right hand side. */
      if (rhs_vals)
        for (i = 0; i < N; i++)
          //FIXME: PUT BACK ORIGINAL:          fprintf (mat_file, "%32.24e\n", rhs_vals[i]);
          fprintf (mat_file, "%32.24e\n", -rhs_vals[i]);

      fclose (mat_file);
    }
  }

  ESymSolverStatus
  IterativePardisoSolverInterface::Factorization(const Index* ia,
      const Index* ja,
      bool check_NegEVals,
      Index numberOfNegEVals)
  {
    DBG_START_METH("IterativePardisoSolverInterface::Factorization",dbg_verbosity);

    // Call Pardiso to do the factorization
    ipfint PHASE ;
    ipfint N = dim_;
    ipfint PERM;   // This should not be accessed by Pardiso
    ipfint NRHS = 0;
    double B;  // This should not be accessed by Pardiso in factorization
    // phase
    double X;  // This should not be accessed by Pardiso in factorization
    // phase
    ipfint ERROR;

    bool done = false;
    bool just_performed_symbolic_factorization = false;

    while (!done) {
      bool is_normal = false;
      if (IsNull(InexData().normal_x()) && InexData().compute_normal()) {
        is_normal = true;
      }
      if (!have_symbolic_factorization_) {
        if (HaveIpData()) {
          IpData().TimingStats().LinearSystemSymbolicFactorization().Start();
        }
        PHASE = 11;
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                       "Calling Pardiso for symbolic factorization.\n");

        if (is_normal) {
          DPARM_[ 0] = normal_pardiso_max_iter_;
          DPARM_[ 1] = normal_pardiso_iter_relative_tol_;
          DPARM_[ 2] = normal_pardiso_iter_coarse_size_;
          DPARM_[ 3] = normal_pardiso_iter_max_levels_;
          DPARM_[ 4] = normal_pardiso_iter_dropping_factor_used_;
          DPARM_[ 5] = normal_pardiso_iter_dropping_schur_used_;
          DPARM_[ 6] = normal_pardiso_iter_max_row_fill_;
          DPARM_[ 7] = normal_pardiso_iter_inverse_norm_factor_;
        }
        else {
          DPARM_[ 0] = pardiso_max_iter_;
          DPARM_[ 1] = pardiso_iter_relative_tol_;
          DPARM_[ 2] = pardiso_iter_coarse_size_;
          DPARM_[ 3] = pardiso_iter_max_levels_;
          DPARM_[ 4] = pardiso_iter_dropping_factor_used_;
          DPARM_[ 5] = pardiso_iter_dropping_schur_used_;
          DPARM_[ 6] = pardiso_iter_max_row_fill_;
          DPARM_[ 7] = pardiso_iter_inverse_norm_factor_;
        }
        DPARM_[ 8] = 25; // maximum number of non-improvement steps

        F77_FUNC(pardiso,PARDISO)(PT_, &MAXFCT_, &MNUM_, &MTYPE_,
                                  &PHASE, &N, a_, ia, ja, &PERM,
                                  &NRHS, IPARM_, &MSGLVL_, &B, &X,
                                  &ERROR, DPARM_);
        if (HaveIpData()) {
          IpData().TimingStats().LinearSystemSymbolicFactorization().End();
        }
        if (ERROR==-7) {
          Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                         "Pardiso symbolic factorization returns ERROR = %d.  Matrix is singular.\n", ERROR);
          return SYMSOLVER_SINGULAR;
        }
        else if (ERROR!=0) {
          Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                         "Error in Pardiso during symbolic factorization phase.  ERROR = %d.\n", ERROR);
          return SYMSOLVER_FATAL_ERROR;
        }
        have_symbolic_factorization_ = true;
        just_performed_symbolic_factorization = true;

        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                       "Memory in KB required for the symbolic factorization  = %d.\n", IPARM_[14]);
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                       "Integer memory in KB required for the numerical factorization  = %d.\n", IPARM_[15]);
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                       "Double  memory in KB required for the numerical factorization  = %d.\n", IPARM_[16]);
      }

      PHASE = 22;

      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemFactorization().Start();
      }
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "Calling Pardiso for factorization.\n");
      // Dump matrix to file, and count number of solution steps.
      if (HaveIpData()) {
        if (IpData().iter_count() != debug_last_iter_)
          debug_cnt_ = 0;
        debug_last_iter_ = IpData().iter_count();
        debug_cnt_ ++;
      }
      else {
        debug_cnt_ = 0;
        debug_last_iter_ = 0;
      }

      if (is_normal) {
        DPARM_[ 0] = normal_pardiso_max_iter_;
        DPARM_[ 1] = normal_pardiso_iter_relative_tol_;
        DPARM_[ 2] = normal_pardiso_iter_coarse_size_;
        DPARM_[ 3] = normal_pardiso_iter_max_levels_;
        DPARM_[ 4] = normal_pardiso_iter_dropping_factor_used_;
        DPARM_[ 5] = normal_pardiso_iter_dropping_schur_used_;
        DPARM_[ 6] = normal_pardiso_iter_max_row_fill_;
        DPARM_[ 7] = normal_pardiso_iter_inverse_norm_factor_;
      }
      else {
        DPARM_[ 0] = pardiso_max_iter_;
        DPARM_[ 1] = pardiso_iter_relative_tol_;
        DPARM_[ 2] = pardiso_iter_coarse_size_;
        DPARM_[ 3] = pardiso_iter_max_levels_;
        DPARM_[ 4] = pardiso_iter_dropping_factor_used_;
        DPARM_[ 5] = pardiso_iter_dropping_schur_used_;
        DPARM_[ 6] = pardiso_iter_max_row_fill_;
        DPARM_[ 7] = pardiso_iter_inverse_norm_factor_;
      }
      DPARM_[ 8] = 25; // maximum number of non-improvement steps

      F77_FUNC(pardiso,PARDISO)(PT_, &MAXFCT_, &MNUM_, &MTYPE_,
                                &PHASE, &N, a_, ia, ja, &PERM,
                                &NRHS, IPARM_, &MSGLVL_, &B, &X,
                                &ERROR, DPARM_);
      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemFactorization().End();
      }

      if (ERROR==-7) {
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                       "Pardiso factorization returns ERROR = %d.  Matrix is singular.\n", ERROR);
        return SYMSOLVER_SINGULAR;
      }
      else if (ERROR==-4) {
        // I think this means that the matrix is singular
        // OLAF said that this will never happen (ToDo)
        return SYMSOLVER_SINGULAR;
      }
      else if (ERROR!=0 ) {
        Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                       "Error in Pardiso during factorization phase.  ERROR = %d.\n", ERROR);
        return SYMSOLVER_FATAL_ERROR;
      }

      negevals_ = Max(IPARM_[22], numberOfNegEVals);
      if (IPARM_[13] != 0) {
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                       "Number of perturbed pivots in factorization phase = %d.\n", IPARM_[13]);
        if (HaveIpData()) {
          IpData().Append_info_string("Pp");
        }
        done = true;
      }
      else {
        done = true;
      }
    }

    //  DBG_ASSERT(IPARM_[21]+IPARM_[22] == dim_);

    // Check whether the number of negative eigenvalues matches the requested
    // count
    if (skip_inertia_check_) numberOfNegEVals=negevals_;

    if (check_NegEVals && (numberOfNegEVals!=negevals_)) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Wrong inertia: required are %d, but we got %d.\n",
                     numberOfNegEVals, negevals_);
      return SYMSOLVER_WRONG_INERTIA;
    }

    return SYMSOLVER_SUCCESS;
  }

  ESymSolverStatus IterativePardisoSolverInterface::Solve(const Index* ia,
      const Index* ja,
      Index nrhs,
      double *rhs_vals)
  {
    DBG_START_METH("IterativePardisoSolverInterface::Solve",dbg_verbosity);

    DBG_ASSERT(nrhs==1);

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemBackSolve().Start();
    }
    // Call Pardiso to do the solve for the given right-hand sides
    ipfint PHASE = 33;
    ipfint N = dim_;
    ipfint PERM;   // This should not be accessed by Pardiso
    ipfint NRHS = nrhs;
    double* X = new double[nrhs*dim_];
    double* ORIG_RHS = new double[nrhs*dim_];
    ipfint ERROR;

    // Initialize solution with zero and save right hand side
    for (int i = 0; i < N; i++) {
      X[i] = 0;
      ORIG_RHS[i] = rhs_vals[i];
    }

    // Dump matrix to file if requested
    Index iter_count = 0;
    if (HaveIpData()) {
      iter_count = IpData().iter_count();
    }
    write_iajaa_matrix (N, ia, ja, a_, rhs_vals, iter_count, debug_cnt_);

    IterativeSolverTerminationTester* tester;

    int attempts = 0;
    const int max_attempts = pardiso_max_droptol_corrections_+1;

    bool is_normal = false;
    if (IsNull(InexData().normal_x()) && InexData().compute_normal()) {
      tester = GetRawPtr(normal_tester_);
      is_normal = true;
    }
    else {
      tester = GetRawPtr(pd_tester_);
    }
    global_tester_ptr_ = tester;

    while (attempts<max_attempts) {
      bool retval = tester->InitializeSolve();
      ASSERT_EXCEPTION(retval, INTERNAL_ABORT, "tester->InitializeSolve(); returned false");

      for (int i = 0; i < N; i++) {
        rhs_vals[i] = ORIG_RHS[i];
      }

      DPARM_[ 8] = 25; // non_improvement in SQMR iteration
      F77_FUNC(pardiso,PARDISO)(PT_, &MAXFCT_, &MNUM_, &MTYPE_,
                                &PHASE, &N, a_, ia, ja, &PERM,
                                &NRHS, IPARM_, &MSGLVL_, rhs_vals, X,
                                &ERROR, DPARM_);

      if (ERROR <= -100 && ERROR >= -110) {
        Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                       "Iterative solver in Pardiso did not converge (ERROR = %d)\n", ERROR);
        Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                       "  Decreasing drop tolerances from DPARM_[ 4] = %e and DPARM_[ 5] = %e ", DPARM_[ 4], DPARM_[ 5]);
        if (is_normal) {
          Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                         "(normal step)\n");
        }
        else {
          Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                         "(PD step)\n");
        }
        PHASE = 23;
        DPARM_[ 4] *= decr_factor_;
        DPARM_[ 5] *= decr_factor_;
        Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                       "                               to DPARM_[ 4] = %e and DPARM_[ 5] = %e\n", DPARM_[ 4], DPARM_[ 5]);
        // Copy solution back to y to get intial values for the next iteration
        attempts++;
        ERROR = 0;
      }
      else  {
        attempts = max_attempts;
        Index iterations_used = tester->GetSolverIterations();
        if (is_normal) {
          Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                         "Number of iterations in Pardiso iterative solver for normal step = %d.\n", iterations_used);
        }
        else {
          Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                         "Number of iterations in Pardiso iterative solver for PD step = %d.\n", iterations_used);
        }
      }
      tester->Clear();
    }

    if (is_normal) {
      if (DPARM_[4] < normal_pardiso_iter_dropping_factor_) {
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                       "Increasing drop tolerances from DPARM_[ 4] = %e and DPARM_[ 5] = %e (normal step\n", DPARM_[ 4], DPARM_[ 5]);
      }
      normal_pardiso_iter_dropping_factor_used_ =
        Min(DPARM_[4]/decr_factor_, normal_pardiso_iter_dropping_factor_);
      normal_pardiso_iter_dropping_schur_used_ =
        Min(DPARM_[5]/decr_factor_, normal_pardiso_iter_dropping_schur_);
      if (DPARM_[4] < normal_pardiso_iter_dropping_factor_) {
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                       "                             to DPARM_[ 4] = %e and DPARM_[ 5] = %e for next iteration.\n", normal_pardiso_iter_dropping_factor_used_, normal_pardiso_iter_dropping_schur_used_);
      }
    }
    else {
      if (DPARM_[4] < pardiso_iter_dropping_factor_) {
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                       "Increasing drop tolerances from DPARM_[ 4] = %e and DPARM_[ 5] = %e (PD step\n", DPARM_[ 4], DPARM_[ 5]);
      }
      pardiso_iter_dropping_factor_used_ =
        Min(DPARM_[4]/decr_factor_, pardiso_iter_dropping_factor_);
      pardiso_iter_dropping_schur_used_ =
        Min(DPARM_[5]/decr_factor_, pardiso_iter_dropping_schur_);
      if (DPARM_[4] < pardiso_iter_dropping_factor_) {
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                       "                             to DPARM_[ 4] = %e and DPARM_[ 5] = %e for next iteration.\n", pardiso_iter_dropping_factor_used_, pardiso_iter_dropping_schur_used_);
      }
    }

    delete [] X;
    delete [] ORIG_RHS;

    if (IPARM_[6] != 0) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Number of iterative refinement steps = %d.\n", IPARM_[6]);
      if (HaveIpData()) {
        IpData().Append_info_string("Pi");
      }
    }

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemBackSolve().End();
    }
    if (ERROR!=0 ) {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "Error in Pardiso during solve phase.  ERROR = %d.\n", ERROR);
      return SYMSOLVER_FATAL_ERROR;
    }
    if (test_result_ == IterativeSolverTerminationTester::MODIFY_HESSIAN) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Termination tester requests modification of Hessian\n");
      return SYMSOLVER_WRONG_INERTIA;
    }
#if 0
    // FRANK: look at this:
    if (test_result_ == IterativeSolverTerminationTester::CONTINUE) {
      if (InexData().compute_normal()) {
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                       "Termination tester not satisfied!!! Pretend singular\n");
        return SYMSOLVER_SINGULAR;
      }
    }
#endif
    if (test_result_ == IterativeSolverTerminationTester::TEST_2_SATISFIED) {
      // Termination Test 2 is satisfied, set the step for the primal
      // iterates to zero
      Index nvars = IpData().curr()->x()->Dim() + IpData().curr()->s()->Dim();
      const Number zero = 0.;
      IpBlasDcopy(nvars, &zero, 0, rhs_vals, 1);
    }
    return SYMSOLVER_SUCCESS;
  }

  Index IterativePardisoSolverInterface::NumberOfNegEVals() const
  {
    DBG_START_METH("IterativePardisoSolverInterface::NumberOfNegEVals",dbg_verbosity);
    DBG_ASSERT(negevals_>=0);
    return negevals_;
  }

  bool IterativePardisoSolverInterface::IncreaseQuality()
  {
    // At the moment, I don't see how we could tell Pardiso to do better
    // (maybe switch from IPARM[20]=1 to IPARM[20]=2?)
    return false;
  }

} // namespace Ipopt
