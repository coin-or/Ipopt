// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter             IBM    2008-09-19
//
//           based on IpPardisoSolverInterface.cpp rev 1219

#include "IpoptConfig.h"
#include "IpIterativePardisoSolverInterface.hpp"
#include "IpBlas.hpp"

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
                                         ipfint* IPARM);
  void F77_FUNC(pardiso,PARDISO)(void** PT, const ipfint* MAXFCT,
                                 const ipfint* MNUM, const ipfint* MTYPE,
                                 const ipfint* PHASE, const ipfint* N,
                                 const double* A, const ipfint* IA,
                                 const ipfint* JA, const ipfint* PERM,
                                 const ipfint* NRHS, ipfint* IPARM,
                                 const ipfint* MSGLVL, double* B, double* X,
                                 ipfint* ERROR);
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
                                &MSGLVL_, &ddmy, &ddmy, &ERROR);
      DBG_ASSERT(ERROR==0);
    }

    delete[] PT_;
    delete[] IPARM_;
    delete[] a_;
  }

  void IterativePardisoSolverInterface::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {}

  bool IterativePardisoSolverInterface::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    Index enum_int;
    options.GetEnumValue("pardiso_matching_strategy", enum_int, prefix);
    match_strat_ = PardisoMatchingStrategy(enum_int);
    options.GetBoolValue("pardiso_redo_symbolic_fact_only_if_inertia_wrong",
                         pardiso_redo_symbolic_fact_only_if_inertia_wrong_,
                         prefix);
    options.GetBoolValue("pardiso_repeated_perturbation_means_singular",
                         pardiso_repeated_perturbation_means_singular_,
                         prefix);
    Index pardiso_out_of_core_power;
    options.GetIntegerValue("pardiso_out_of_core_power",
                            pardiso_out_of_core_power, prefix);
    options.GetBoolValue("pardiso_skip_inertia_check",
                         skip_inertia_check_, prefix);
    options.GetIntegerValue("pardiso_dropping_schur_exponent",
                            pardiso_dropping_schur_exponent_, prefix);
    options.GetIntegerValue("pardiso_dropping_factor_exponent",
                            pardiso_dropping_factor_exponent_, prefix);
    options.GetIntegerValue("pardiso_inverse_norm_factor",
                            pardiso_inverse_norm_factor_, prefix);
    options.GetIntegerValue("pardiso_max_iter", pardiso_max_iter_, prefix);
    options.GetIntegerValue("pardiso_dropping_schur_exponent",
                            normal_pardiso_dropping_schur_exponent_,
                            prefix+"normal.");
    options.GetIntegerValue("pardiso_dropping_factor_exponent",
                            normal_pardiso_dropping_factor_exponent_,
                            prefix+"normal.");
    options.GetIntegerValue("pardiso_inverse_norm_factor",
                            normal_pardiso_inverse_norm_factor_,
                            prefix+"normal.");
    options.GetIntegerValue("pardiso_max_iter", normal_pardiso_max_iter_,
                            prefix+"normal.");
    int pardiso_msglvl;
    options.GetIntegerValue("pardiso_msglvl", pardiso_msglvl, prefix);

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
                                &MSGLVL_, &ddmy, &ddmy, &ERROR);
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
    F77_FUNC(pardisoinit,PARDISOINIT)(PT_, &MTYPE_, IPARM_);

    // Set some parameters for Pardiso
    IPARM_[0] = 1;  // Don't use the default values

#ifdef HAVE_PARDISO_PARALLEL
    // Obtain the numbers of processors from the value of OMP_NUM_THREADS
    char    *var = getenv("OMP_NUM_THREADS");
    int      num_procs;
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
    else {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "You need to set environment variable OMP_NUM_THREADS to the number of processors used in Pardiso (e.g., 1).\n\n");
      return false;
    }
    IPARM_[2] = num_procs;  // Set the number of processors
#else

    IPARM_[2] = 1;
#endif

    IPARM_[1] = 5;
    IPARM_[5] = 1;  // Overwrite right-hand side
    // ToDo: decide if we need iterative refinement in Pardiso.  For
    // now, switch it off ?
    IPARM_[7] = 0;

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
    IPARM_[29] = 1; // we need this for IPOPT interface

    IPARM_[39] = 10 ;  // it was 4 max fill for factor
    IPARM_[40] = 1 ;  // mantisse dropping value for schur complement
    IPARM_[46] = 1 ;  // mantisse dropping value for incomplete factor
    IPARM_[48] = 1 ;  // iterative solver
    MSGLVL_ = pardiso_msglvl;

    // Option for the out of core variant
    IPARM_[49] = pardiso_out_of_core_power;

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

      sprintf (mat_name, "%s_%03d-%02d.iajaa", mat_pref, iter_cnt, sol_cnt);

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
          fprintf (mat_file, "%32.24e\n", rhs_vals[i]);

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
      if (!have_symbolic_factorization_) {
        if (HaveIpData()) {
          IpData().TimingStats().LinearSystemSymbolicFactorization().Start();
        }
        PHASE = 11;
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                       "Calling Pardiso for symbolic factorization.\n");
        F77_FUNC(pardiso,PARDISO)(PT_, &MAXFCT_, &MNUM_, &MTYPE_,
                                  &PHASE, &N, a_, ia, ja, &PERM,
                                  &NRHS, IPARM_, &MSGLVL_, &B, &X,
                                  &ERROR);
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

      bool is_normal = false;
      if (!IsValid(InexData().normal_x())) {
        is_normal = true;
      }
      if (is_normal) {
        IPARM_[41] = normal_pardiso_dropping_schur_exponent_;
        // it  exponent dropping value for schur complement
        IPARM_[42] = normal_pardiso_max_iter_; // max number of iterations
        IPARM_[43] = normal_pardiso_inverse_norm_factor_; // norm of the inverse for algebraic solver
        IPARM_[44] = normal_pardiso_dropping_factor_exponent_;  // exponent dropping value for incomplete factor
      }
      else {
        IPARM_[41] = pardiso_dropping_schur_exponent_;
        // it  exponent dropping value for schur complement
        IPARM_[42] = pardiso_max_iter_; // max number of iterations
        IPARM_[43] = pardiso_inverse_norm_factor_; // norm of the inverse for algebraic solver
        IPARM_[44] = pardiso_dropping_factor_exponent_;  // exponent dropping value for incomplete factor
      }

      F77_FUNC(pardiso,PARDISO)(PT_, &MAXFCT_, &MNUM_, &MTYPE_,
                                &PHASE, &N, a_, ia, ja, &PERM,
                                &NRHS, IPARM_, &MSGLVL_, &B, &X,
                                &ERROR);
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
    ipfint ERROR;

    // Dump matrix to file if requested
    Index iter_count = 0;
    if (HaveIpData()) {
      iter_count = IpData().iter_count();
    }
    write_iajaa_matrix (N, ia, ja, a_, rhs_vals, iter_count, debug_cnt_);

    IterativeSolverTerminationTester* tester;
    bool is_normal = false;
    if (IsNull(InexData().normal_x()) && InexData().compute_normal()) {
      tester = GetRawPtr(normal_tester_);
      is_normal = true;
    }
    else {
      tester = GetRawPtr(pd_tester_);
    }
    global_tester_ptr_ = tester;

    bool retval = tester->InitializeSolve();
    ASSERT_EXCEPTION(retval, INTERNAL_ABORT, "tester->InitializeSolve(); returned false");

    F77_FUNC(pardiso,PARDISO)(PT_, &MAXFCT_, &MNUM_, &MTYPE_,
                              &PHASE, &N, a_, ia, ja, &PERM,
                              &NRHS, IPARM_, &MSGLVL_, rhs_vals, X,
                              &ERROR);

    delete [] X; /* OLAF/MICHAEL: do we really need X? */

    Index iterations_used = tester->GetSolverIterations();
    if (is_normal) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Number of iterations in Pardiso iterative solver for normal step = %d.\n", iterations_used);
    }
    else {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Number of iterations in Pardiso iterative solver for PD step = %d.\n", iterations_used);
    }

    tester->Clear();

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
    // FRANK: look at this:
    if (test_result_ == IterativeSolverTerminationTester::CONTINUE) {
      if (InexData().compute_normal()) {
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                       "Termination tester not satisfied!!! Pretend singular\n");
        return SYMSOLVER_SINGULAR;
      }/*
            else {
       Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
               "Termination tester not satisfied!!! Compute normal step\n");
              InexData().set_next_compute_normal(true);
       return SYMSOLVER_SUCCESS; // Will setting success work out in the end?
       }*/
    }
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
