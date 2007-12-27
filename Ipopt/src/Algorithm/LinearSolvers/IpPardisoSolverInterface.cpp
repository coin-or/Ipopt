// Copyright (C) 2005, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-03-17
//
//           Olaf Schenk                      Univ of Basel 2005-09-20
//                  - changed options, added PHASE_ flag

#include "IpoptConfig.h"
#include "IpPardisoSolverInterface.hpp"

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

  PardisoSolverInterface::PardisoSolverInterface()
      :
      a_(NULL),
      negevals_(-1),
      initialized_(false),

      MAXFCT_(1),
      MNUM_(1),
      MTYPE_(-2),
      MSGLVL_(0),
      debug_last_iter_(-1)
  {
    DBG_START_METH("PardisoSolverInterface::PardisoSolverInterface()",dbg_verbosity);

    PT_ = new void*[64];
    IPARM_ = new ipfint[64];
  }

  PardisoSolverInterface::~PardisoSolverInterface()
  {
    DBG_START_METH("PardisoSolverInterface::~PardisoSolverInterface()",
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

  void PardisoSolverInterface::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    // Todo Use keywords instead of integer numbers
    roptions->AddStringOption3(
      "pardiso_matching_strategy",
      "Matching strategy to be used by Pardiso",
      "complete+2x2",
      "complete", "Match complete (IPAR(13)=1)",
      "complete+2x2", "Match complete+2x2 (IPAR(13)=2)",
      "constraints", "Match constraints (IPAR(13)=3)",
      "This is IPAR(13) in Pardiso manual.  This option is only available if "
      "Ipopt has been compiled with Pardiso.");
    roptions->AddStringOption2(
      "pardiso_redo_symbolic_fact_only_if_inertia_wrong",
      "Toggel for handling case when elements were pertured by Pardiso.",
      "no",
      "no", "Always redo symbolic factorization when elements were perturbed",
      "yes", "Only redo symbolic factorization when elements were perturbed if also the inertia was wrong",
      "This option is only available if Ipopt has been compiled with Pardiso.");
    roptions->AddStringOption2(
      "pardiso_repeated_perturbation_means_singular",
      "Interpretation of perturbed elements.",
      "no",
      "no", "Don't assume that matrix is singular if elements were perturbed after recent symbolic factorization",
      "yes", "Assume that matrix is singular if elements were perturbed after recent symbolic factorization",
      "This option is only available if Ipopt has been compiled with Pardiso.");
    roptions->AddLowerBoundedIntegerOption(
      "pardiso_out_of_core_power",
      "Enables out-of-core variant of Pardiso",
      0, 0,
      "Setting this option to a positive integer k makes Pardiso work in the "
      "out-of-core variant where the factor is split in 2^k subdomains.  This "
      "is IPARM(50) in the Pardiso manual.  This option is only available if "
      "Ipopt has been compiled with Pardiso.");
    roptions->AddStringOption2(
      "pardiso_skip_inertia_check",
      "Always pretent inertia is correct.",
      "no",
      "no", "check interia",
      "yes", "skip inertia check",
      "Setting this option to \"yes\" essentially disables inertia check. "
      "This option makes the algorithm non-robust and easily fail, but it "
      "might give some insight into the necessity of interia control.");
    roptions->AddIntegerOption(
      "pardiso_iter_tol_exponent",
      "",
      -14,
      "");
    roptions->AddStringOption2(
      "pardiso_iterative",
      "",
      "no",
      "no", "",
      "yes", "",
      "");
  }

  bool PardisoSolverInterface::InitializeImpl(const OptionsList& options,
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
    bool pardiso_iterative;
    options.GetBoolValue("pardiso_iterative", pardiso_iterative,
                         prefix);
    int pardiso_iter_tol_exponent;
    options.GetIntegerValue("pardiso_iter_tol_exponent",
                            pardiso_iter_tol_exponent, prefix);

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

    IPARM_[39] = 4 ;  // it was 4 max fill for factor
    IPARM_[40] = 1 ;  // mantisse dropping value for schur complement
    IPARM_[41] =-3 ;  // it  exponent dropping value for schur complement
    IPARM_[42] = 200; // max number of iterations
    IPARM_[43] = 500 ; // norm of the inverse for algebraic solver
    IPARM_[44] =-3 ;  // exponent dropping value for incomplete factor
    IPARM_[46] = 1 ;  // mantisse dropping value for incomplete factor
    IPARM_[45] = pardiso_iter_tol_exponent ;  // residual tolerance
    IPARM_[48] = pardiso_iterative ? 1 : 0 ;  // active direct solver
    if (pardiso_iterative) MSGLVL_ = 2;

    // Option for the out of core variant
    IPARM_[49] = pardiso_out_of_core_power;

    return true;
  }

  ESymSolverStatus PardisoSolverInterface::MultiSolve(bool new_matrix,
      const Index* ia,
      const Index* ja,
      Index nrhs,
      double* rhs_vals,
      bool check_NegEVals,
      Index numberOfNegEVals)
  {
    DBG_START_METH("PardisoSolverInterface::MultiSolve",dbg_verbosity);
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

  double* PardisoSolverInterface::GetValuesArrayPtr()
  {
    DBG_ASSERT(initialized_);
    DBG_ASSERT(a_);
    return a_;
  }

  /** Initialize the local copy of the positions of the nonzero
      elements */
  ESymSolverStatus PardisoSolverInterface::InitializeStructure
  (Index dim, Index nonzeros,
   const Index* ia,
   const Index* ja)
  {
    DBG_START_METH("PardisoSolverInterface::InitializeStructure",dbg_verbosity);
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
  PardisoSolverInterface::SymbolicFactorization(const Index* ia,
      const Index* ja)
  {
    DBG_START_METH("PardisoSolverInterface::SymbolicFactorization",
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
  PardisoSolverInterface::Factorization(const Index* ia,
                                        const Index* ja,
                                        bool check_NegEVals,
                                        Index numberOfNegEVals)
  {
    DBG_START_METH("PardisoSolverInterface::Factorization",dbg_verbosity);

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
        if ( !pardiso_redo_symbolic_fact_only_if_inertia_wrong_ ||
             (negevals_ != numberOfNegEVals) ) {
          if (HaveIpData()) {
            IpData().Append_info_string("Pn");
          }
          have_symbolic_factorization_ = false;
          // We assume now that if there was just a symbolic
          // factorization and we still have perturbed pivots, that
          // the system is actually singular, if
          // pardiso_repeated_perturbation_means_singular_ is true
          if (just_performed_symbolic_factorization) {
            if (pardiso_repeated_perturbation_means_singular_) {
              if (HaveIpData()) {
                IpData().Append_info_string("Ps");
              }
              return SYMSOLVER_SINGULAR;
            }
            else {
              done = true;
            }
          }
          else {
            done = false;
          }
        }
        else {
          if (HaveIpData()) {
            IpData().Append_info_string("Pp");
          }
          done = true;
        }
      }
      else {
        done = true;
      }
    }

    DBG_ASSERT(IPARM_[21]+IPARM_[22] == dim_);

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

  ESymSolverStatus PardisoSolverInterface::Solve(const Index* ia,
      const Index* ja,
      Index nrhs,
      double *rhs_vals)
  {
    DBG_START_METH("PardisoSolverInterface::Solve",dbg_verbosity);

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

    F77_FUNC(pardiso,PARDISO)(PT_, &MAXFCT_, &MNUM_, &MTYPE_,
                              &PHASE, &N, a_, ia, ja, &PERM,
                              &NRHS, IPARM_, &MSGLVL_, rhs_vals, X,
                              &ERROR);

    delete [] X; /* OLAF/MICHAEL: do we really need X? */

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
    return SYMSOLVER_SUCCESS;
  }

  Index PardisoSolverInterface::NumberOfNegEVals() const
  {
    DBG_START_METH("PardisoSolverInterface::NumberOfNegEVals",dbg_verbosity);
    DBG_ASSERT(negevals_>=0);
    return negevals_;
  }

  bool PardisoSolverInterface::IncreaseQuality()
  {
    // At the moment, I don't see how we could tell Pardiso to do better
    // (maybe switch from IPARM[20]=1 to IPARM[20]=2?)
    return false;
  }

} // namespace Ipopt
