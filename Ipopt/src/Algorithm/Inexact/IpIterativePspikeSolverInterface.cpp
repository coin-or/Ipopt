// Copyright (C) 2008, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Madan Sathe (University of Basel), Andreas Waechter (IBM)   2009-08-17
//
//           based on IpIterativePardisoSolverInterface.cpp rev 1535

#include "IpoptConfig.h"
#include "IpIterativePspikeSolverInterface.hpp"
#include "IpBlas.hpp"
# include <math.h>

#include "IpMpi.hpp"

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

extern Ipopt::IterativeSolverTerminationTester* global_tester_ptr_;
extern Ipopt::IterativeSolverTerminationTester::ETerminationTest test_result_;
extern "C"
{
  int IpoptTerminationTest(int n, double* sol, double* resid, int iter, double norm2_rhs);

  // The following global function pointer is defined in the Pardiso library
  void SetIpoptCallbackFunction(int (*IpoptFunction)(int n, double* x,  double* r, int k, double b));
}

/** Prototypes for Pardiso's subroutines */
extern "C"
{
  // INIT PSPIKE: job = 0, FACT PSPIKE = 1, SOLVE PSPIKE = 2, DESTROY PSPIKE = 3
  void F77_FUNC(pspike, PSPIKE)(int * job, int * neqns, int * nzmax,
                                const ipfint * ia, const ipfint * ja, const double * a,
                                double * f, int * kk, double * tol, int * nrhs);
}

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  IterativePspikeSolverInterface::
  IterativePspikeSolverInterface(IterativeSolverTerminationTester& normal_tester,
                                 IterativeSolverTerminationTester& pd_tester)
      :
      a_(NULL),
      initialized_(false),

      debug_last_iter_(-1),
      normal_tester_(&normal_tester),
      pd_tester_(&pd_tester)
  {
    DBG_START_METH("IterativePspikeSolverInterface::IterativePspikeSolverInterface()",dbg_verbosity);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_);

    int neqns = 0;
    int nzmax = 0;
    const ipfint * ia = NULL;
    const ipfint * ja = NULL;
    const double * a = NULL;
    double * f = NULL;
    int bandwidth = 0;
    double tol = 0.0;
    int nrhs = 0;

    int job = 0;
    F77_FUNC(pspike,PSPIKE)(&job, &neqns, &nzmax, ia, ja, a, f, &bandwidth, &tol, &nrhs);
  }

  IterativePspikeSolverInterface::~IterativePspikeSolverInterface()
  {
    DBG_START_METH("IterativePspikeSolverInterface::~IterativePspikeSolverInterface()",
                   dbg_verbosity);
    delete[] a_;

    if (initialized_) {
      int neqns = 0;
      int nzmax = 0;
      const ipfint * ia = NULL;
      const ipfint * ja = NULL;
      const double * a = NULL;
      double * f = NULL;
      int bandwidth = 0;
      double tol = 0.0;
      int nrhs = 0;

      int job = 4;
      F77_FUNC(pspike,PSPIKE)(&job, &neqns, &nzmax, ia, ja, a, f, &bandwidth, &tol, &nrhs);
      // Madan: Do we need to initialize again afterwards?
    }
  }

  void IterativePspikeSolverInterface::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddLowerBoundedNumberOption(
      "pspike_tol",
      "",
      0.0, true,
      1e-4,
      "");
    roptions->AddLowerBoundedIntegerOption(
      "pspike_bandwidth",
      "",
      0, 100,
      "");
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


  bool IterativePspikeSolverInterface::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    options.GetNumericValue("pspike_tol", pspike_tol_, prefix);
    options.GetIntegerValue("pspike_bandwidth", pspike_bandwidth_, prefix);

    // DESTROY PSPIKE with job = 4 !!!!
    if (initialized_) {
      int neqns = 0;
      int nzmax = 0;
      const ipfint * ia = NULL;
      const ipfint * ja = NULL;
      const double * a = NULL;
      double * f = NULL;
      int bandwidth = 0;
      double tol = 0.0;
      int nrhs = 0;

      int job = 4;
      F77_FUNC(pspike,PSPIKE)(&job, &neqns, &nzmax, ia, ja, a, f, &bandwidth, &tol, &nrhs);
      job = 0;
      F77_FUNC(pspike,PSPIKE)(&job, &neqns, &nzmax, ia, ja, a, f, &bandwidth, &tol, &nrhs);
      // Madan: Do we need to initialize again afterwards?
    }

    have_symbolic_factorization_ = false;

    SetIpoptCallbackFunction(&IpoptTerminationTest);

    bool retval = normal_tester_->Initialize(Jnlst(), IpNLP(), IpData(),
                  IpCq(), options, prefix);
    if (retval) {
      retval = pd_tester_->Initialize(Jnlst(), IpNLP(), IpData(),
                                      IpCq(), options, prefix);
    }

    return retval;
  }

  ESymSolverStatus IterativePspikeSolverInterface::MultiSolve(bool new_matrix,
      const Index* ia,
      const Index* ja,
      Index nrhs,
      double* rhs_vals,
      bool check_NegEVals,
      Index numberOfNegEVals)
  {
    DBG_START_METH("IterativePspikeSolverInterface::MultiSolve",dbg_verbosity);
    DBG_ASSERT(!check_NegEVals || ProvidesInertia());
    DBG_ASSERT(initialized_);

    if (my_rank_ == 0) {
      write_iajaa_matrix (dim_, ia, ja, a_, rhs_vals, 0, 0);
    }

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
    ESymSolverStatus status = Solve(ia, ja, nrhs, rhs_vals);

    return status;
  }

  double* IterativePspikeSolverInterface::GetValuesArrayPtr()
  {
    DBG_ASSERT(initialized_);
    DBG_ASSERT(a_);
    return a_;
  }

  /** Initialize the local copy of the positions of the nonzero
      elements */
  ESymSolverStatus IterativePspikeSolverInterface::InitializeStructure
  (Index dim, Index nonzeros,
   const Index* ia,
   const Index* ja)
  {
    DBG_START_METH("IterativePspikeSolverInterface::InitializeStructure",dbg_verbosity);
    dim_ = dim;
    nonzeros_ = nonzeros;

    if (my_rank_==0) {
      // Make space for storing the matrix elements
      delete[] a_;
      a_ = NULL;
      a_ = new double[nonzeros_];
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
  IterativePspikeSolverInterface::SymbolicFactorization(const Index* ia,
      const Index* ja)
  {
    DBG_START_METH("IterativePspikeSolverInterface::SymbolicFactorization",
                   dbg_verbosity);

    return SYMSOLVER_SUCCESS;
  }



  ESymSolverStatus
  IterativePspikeSolverInterface::Factorization(const Index* ia,
      const Index* ja,
      bool check_NegEVals,
      Index numberOfNegEVals)
  {
    DBG_START_METH("IterativePspikeSolverInterface::Factorization",dbg_verbosity);

    double tol = pspike_tol_;
    int bandwidth = pspike_bandwidth_;
    int nzmax = 0;

    if (my_rank_==0) {
      nzmax = ia[dim_]-1;
    }

    if (!have_symbolic_factorization_) {
      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemSymbolicFactorization().Start();
      }

      // ORDERING in PSPIKE with job = 1 !!!!
      int job = 1;

      double * rhs_vals = NULL;
      int nrhs = 1;
      F77_FUNC(pspike,PSPIKE)(&job, &dim_, &nzmax, ia, ja, a_, rhs_vals, &bandwidth, &tol, &nrhs);

      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemSymbolicFactorization().End();
      }
      have_symbolic_factorization_ = true;
    }

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemFactorization().Start();
    }

    // FACTORIZE in PSPIKE with job = 2 !!!!
    int job = 2;

    double * rhs_vals = NULL;
    int nrhs = 1;
    F77_FUNC(pspike,PSPIKE)(&job, &dim_, &nzmax, ia, ja, a_, rhs_vals, &bandwidth, &tol, &nrhs);

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemFactorization().End();
    }

    return SYMSOLVER_SUCCESS;
  }

  ESymSolverStatus IterativePspikeSolverInterface::Solve(const Index* ia,
      const Index* ja,
      Index nrhs,
      double *rhs_vals)
  {
    DBG_START_METH("IterativePspikeSolverInterface::Solve",dbg_verbosity);

    DBG_ASSERT(nrhs==1);

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemBackSolve().Start();
    }

    ipfint N = dim_;
    ipfint NRHS = nrhs;
    ipfint ERROR;
    double* ORIG_RHS = NULL;

    if (my_rank_ == 0) {
      ORIG_RHS = new double[nrhs*dim_];

      // Initialize solution with zero and save right hand side
      for (int i = 0; i < N; i++) {
        ORIG_RHS[i] = rhs_vals[i];
      }

      // Dump matrix to file if requested
      Index iter_count = 0;
      if (HaveIpData()) {
        iter_count = IpData().iter_count();
      }

      write_iajaa_matrix (N, ia, ja, a_, rhs_vals, iter_count, debug_cnt_);
    }

    IterativeSolverTerminationTester* tester;

    int attempts = 0;
    const int max_attempts = 100; //TODO: pardiso_max_droptol_corrections_+1;

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

      double tol = pspike_tol_;
      int bandwidth = pspike_bandwidth_;
      int nzmax = -1;

      if (my_rank_ == 0) {
        for (int i = 0; i < N; i++) {
          rhs_vals[i] = ORIG_RHS[i];
        }
        nzmax = ia[dim_]-1;
      }

      // SOLVE in PSPIKE with job = 3
      int job = 3;
      F77_FUNC(pspike,PSPIKE)(&job, &N, &nzmax, ia, ja, a_, rhs_vals, &bandwidth, &tol, &NRHS);
      ERROR = 0;
      attempts = max_attempts;
      Index iterations_used = tester->GetSolverIterations();
    }
    tester->Clear();

    delete [] ORIG_RHS;

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemBackSolve().End();
    }
    if (ERROR!=0 ) {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "Error in Pspike during solve phase.  ERROR = %d.\n", ERROR);
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
      }
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

  Index IterativePspikeSolverInterface::NumberOfNegEVals() const
  {
    DBG_START_METH("IterativePspikeSolverInterface::NumberOfNegEVals",dbg_verbosity);
    DBG_ASSERT(false && "We should not get here");
    return -1;
  }

  bool IterativePspikeSolverInterface::IncreaseQuality()
  {
    // At the moment, I don't see how we could tell Pardiso to do better
    // (maybe switch from IPARM[20]=1 to IPARM[20]=2?)
    return false;
  }

} // namespace Ipopt
