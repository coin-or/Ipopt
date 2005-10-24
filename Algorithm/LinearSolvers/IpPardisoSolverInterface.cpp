// Copyright (C) 2005 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-03-17
//
//           Olaf Schenk                      Univ of Basel 2005-09-20
//                  - changed options, added PHASE_ flag

#include "IpPardisoSolverInterface.hpp"


#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
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
#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  PardisoSolverInterface::PardisoSolverInterface()
      :
      a_(NULL),
      initialized_(false),

      MAXFCT_(1),
      MNUM_(1),
      MTYPE_(-2),
      MSGLVL_(0)
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
      F77_FUNC(pardiso,PARDISO)(PT_, &MAXFCT_, &MNUM_, &MTYPE_, &PHASE, &N,
                                NULL, NULL, NULL, NULL, &NRHS, IPARM_,
                                &MSGLVL_, NULL, NULL, &ERROR);
      DBG_ASSERT(ERROR==0);
    }

    delete[] PT_;
    delete[] IPARM_;
    delete[] a_;
  }

  bool PardisoSolverInterface::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    // Number value = 0.0;

    // Tell Pardiso to release all memory if it had been used before
    if (initialized_) {
      ipfint PHASE = -1;
      ipfint N = dim_;
      ipfint NRHS = 0;
      ipfint ERROR;
      F77_FUNC(pardiso,PARDISO)(PT_, &MAXFCT_, &MNUM_, &MTYPE_, &PHASE, &N,
                                NULL, NULL, NULL, NULL, &NRHS, IPARM_,
                                &MSGLVL_, NULL, NULL, &ERROR);
      DBG_ASSERT(ERROR==0);
    }

    // Reset all private data
    dim_=0;
    nonzeros_=0;
    have_symbolic_factorization_=false;
    initialized_=false;
    delete[] a_;

    // Call Pardiso's initialization routine
    IPARM_[0] = 0;  // Tell it to fill IPARM with default values(?)
    F77_FUNC(pardisoinit,PARDISOINIT)(PT_, &MTYPE_, IPARM_);

    // Set some parameters for Pardiso
    IPARM_[0] = 1;  // Don't use the default values
    IPARM_[2] = 1;  // Only one CPU for now
    IPARM_[5] = 1;  // Overwrite right-hand side
    // ToDo: decide if we need iterative refinement in Pardiso.  For
    // now, switch it off ?
    IPARM_[7] = 0;

    // Options suggested by Olaf Schenk
    IPARM_[9] = 12;
    IPARM_[10] = 1;
    // Matching information:  IPARM_[12] = 1 seems ok, but results in a
    // large number of pivot perturbation
    // Matching information:  IPARM_[12] = 2 robust,  but more  expensive method
    IPARM_[12] = 2;

    IPARM_[20] = 1;

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

  ESymSolverStatus
  PardisoSolverInterface::Factorization(const Index* ia,
                                        const Index* ja,
                                        bool check_NegEVals,
                                        Index numberOfNegEVals)
  {
    DBG_START_METH("PardisoSolverInterface::Factorization",dbg_verbosity);

    // Call Pardiso to do the factorization
    ipfint N = dim_;
    ipfint PHASE;
    ipfint PERM;   // This should not be accessed by Pardiso
    ipfint NRHS = 0;
    double B;  // This should not be accessed by Pardiso in factorization
    // phase
    double X;  // This should not be accessed by Pardiso in factorization
    // phase
    ipfint ERROR;

    if (getenv ("IPOPT_WRITE_MAT")) {

      // Dump matrix to file...
      ipfint  NNZ = ia[N]-1;
      ipfint   i;

      if (IpData().iter_count() != debug_last_iter_) {
        debug_cnt_ = 0;
      }

      /* Write header */
      FILE    *mat_file;
      char     mat_name[128];
      char     mat_pref[32];

      if (getenv ("IPOPT_WRITE_PREFIX"))
        strcpy (mat_pref, getenv ("IPOPT_WRITE_PREFIX"));
      else
        strcpy (mat_pref, "mat-ipopt");

      sprintf (mat_name, "%s_%03d-%02d.iajaa", mat_pref, IpData().iter_count(), debug_cnt_);
      mat_file = fopen (mat_name, "w");

      fprintf (mat_file, "%d\n", N);
      fprintf (mat_file, "%d\n", NNZ);

      /* Write ia's */
      for (i = 0; i < N+1; i++)
        fprintf (mat_file, "%d\n", ia[i]);

      /* Write ja's */
      for (i = 0; i < NNZ; i++)
        fprintf (mat_file, "%d\n", ja[i]);

      /* Write values */
      for (i = 0; i < NNZ; i++)
        fprintf (mat_file, "%32.24e\n", a_[i]);

      // /* Write RHS */
      // for (i = 0; i < N; i++)
      // 	fprintf (mat_file, "%32.24e\n", rhs_vals[i]);

      fclose (mat_file);

      debug_last_iter_ = IpData().iter_count();
      debug_cnt_++;
    }

    bool done = false;
    bool just_performed_symbolic_factorization = false;

    while (!done) {
      if (!have_symbolic_factorization_) {
        IpData().TimingStats().LinearSystemSymbolicFactorization().Start();
        PHASE = 11;
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                       "Calling Pardiso for symbolic factorization.\n");
        F77_FUNC(pardiso,PARDISO)(PT_, &MAXFCT_, &MNUM_, &MTYPE_,
                                  &PHASE, &N, a_, ia, ja, &PERM,
                                  &NRHS, IPARM_, &MSGLVL_, &B, &X,
                                  &ERROR);
        IpData().TimingStats().LinearSystemSymbolicFactorization().End();
        if (ERROR!=0 ) {
          Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                         "Error in Pardiso during symbolic factorization phase.  ERROR = %d.\n", ERROR);
          return SYMSOLVER_FATAL_ERROR;
        }
        have_symbolic_factorization_ = true;
        just_performed_symbolic_factorization = true;
      }

      PHASE = 22;

      IpData().TimingStats().LinearSystemFactorization().Start();
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "Calling Pardiso for factorization.\n");
      F77_FUNC(pardiso,PARDISO)(PT_, &MAXFCT_, &MNUM_, &MTYPE_,
                                &PHASE, &N, a_, ia, ja, &PERM,
                                &NRHS, IPARM_, &MSGLVL_, &B, &X,
                                &ERROR);
      IpData().TimingStats().LinearSystemFactorization().End();

      if (ERROR==-4) {
        // I think this means that the matrix is singular
        // OLAF said that this will never happen (ToDo)
        return SYMSOLVER_SINGULAR;
      }
      else if (ERROR!=0 ) {
        Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                       "Error in Pardiso during factorization phase.  ERROR = %d.\n", ERROR);
        return SYMSOLVER_FATAL_ERROR;
      }

      negevals_ = IPARM_[22];
      if (IPARM_[13] != 0) {
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                       "Number of perturbed pivots in factorization phase = %d.\n", IPARM_[13]);
        // trigger a new symblic factorization for the next factorize
        // call, even if the current inertia estimate is correct
        // OLAF MICHAEL : Maybe we need to discuss this
        have_symbolic_factorization_ = false;
        // We assume now that if there was just a symbolic
        // factorization and we still have perturbed pivots, that the
        // system is actually singular
        if (just_performed_symbolic_factorization) {
          IpData().Append_info_string("Ps");
          //break; // Declaring this system singular doesn't seem to be
          // working well???
          return SYMSOLVER_SINGULAR;
        }
        IpData().Append_info_string("Pp");
        /*
               // If we there were perturbed pivots, and the inertia of the
               // system is correct, and we haven't redone the symblic
               // factorization during this call yet, trigger a new
               // factorization after a symblic factorization
               done = (just_performed_symbolic_factorization || !check_NegEVals ||
                       numberOfNegEVals==negevals_);
        */
        done = false;
      }
      else {
        done = true;
      }
    }

    DBG_ASSERT(IPARM_[21]+IPARM_[22] == dim_);

    // Check whether the number of negative eigenvalues matches the requested
    // count
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

    IpData().TimingStats().LinearSystemBackSolve().Start();
    // Call Pardiso to do the solve for the given right-hand sides
    ipfint PHASE = 33;
    ipfint N = dim_;
    ipfint PERM;   // This should not be accessed by Pardiso
    ipfint NRHS = nrhs;
    double* X = new double[nrhs*dim_];
    ipfint ERROR;

    F77_FUNC(pardiso,PARDISO)(PT_, &MAXFCT_, &MNUM_, &MTYPE_,
                              &PHASE, &N, a_, ia, ja, &PERM,
                              &NRHS, IPARM_, &MSGLVL_, rhs_vals, X,
                              &ERROR);

    // The following delete was missing (memory leak?)
    delete [] X; /* OLAF/MICHAEL: do we really need X? */

    if (IPARM_[6] != 0) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Number of iterative refinement steps = %d.\n", IPARM_[6]);
      IpData().Append_info_string("Pi");
    }

    IpData().TimingStats().LinearSystemBackSolve().End();
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
