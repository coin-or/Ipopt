// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpPardisoSolverInterface.hpp"

#ifdef OLD_C_HEADERS
# include <math.h>
#else
# include <cmath>
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

  DBG_SET_VERBOSITY(0);

  PardisoSolverInterface::PardisoSolverInterface()
      :
      dim_(0),
      nonzeros_(0),
      initialized_(false),
      negevals_(-1),

      MAXFCT_(1),
      MNUM_(1),
      MTYPE_(-2),
      MSGLVL_(0),

      a_(NULL)
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
    Number value = 0.0;

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

    // IPARM_[20] = 2;

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
    if (retval != SYMSOLVER_SUCCESS ) {
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

    // Call Pardiso to do the analysis phase
    ipfint PHASE = 11;
    ipfint N = dim_;
    ipfint PERM;   // This should not be accessed by Pardiso
    ipfint NRHS = 0;
    double B;  // This should not be accessed by Pardiso in analysis
	       // phase
    double X;  // This should not be accessed by Pardiso in analysis
	       // phase
    ipfint ERROR;

    F77_FUNC(pardiso,PARDISO)(PT_, &MAXFCT_, &MNUM_, &MTYPE_,
			      &PHASE, &N, a_, ia, ja, &PERM,
			      &NRHS, IPARM_, &MSGLVL_, &B, &X,
			      &ERROR);
    if (ERROR!=0) {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "Error in Pardiso during analysis phase.  ERROR = %d.\n",
		     ERROR);
      return SYMSOLVER_FATAL_ERROR;
    }

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
    ipfint PHASE = 22;
    ipfint N = dim_;
    ipfint PERM;   // This should not be accessed by Pardiso
    ipfint NRHS = 0;
    double B;  // This should not be accessed by Pardiso in factorization
	       // phase
    double X;  // This should not be accessed by Pardiso in factorization
	       // phase
    ipfint ERROR;

    F77_FUNC(pardiso,PARDISO)(PT_, &MAXFCT_, &MNUM_, &MTYPE_,
			      &PHASE, &N, a_, ia, ja, &PERM,
			      &NRHS, IPARM_, &MSGLVL_, &B, &X,
			      &ERROR);
    if (ERROR==-4) {
      // I think this means that the matrix is singular
      return SYMSOLVER_SINGULAR;
    }
    else if (ERROR!=0 ) {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "Error in Pardiso during factorization phase.  ERROR = %d.\n", ERROR);
      return SYMSOLVER_FATAL_ERROR;
    }

    /* ToDo ask Olaf what this means
    if (IPARM_[13] != 0) {
      Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                     "Number of perturbed pivots in factorization phase = %d.\n", IPARM_[13]);
    }
    */

    negevals_ = IPARM_[22];

    DBG_ASSERT(IPARM_[21]+IPARM_[22] == dim_);

    // Check whether the number of negative eigenvalues matches the requested
    // count
    if (check_NegEVals && (numberOfNegEVals!=negevals_)) {
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
    if (ERROR!=0 ) {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "Error in Pardiso during factorization phase.  ERROR = %d.\n", ERROR);
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
