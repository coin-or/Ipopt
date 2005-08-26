// Copyright (C) 2005 Yifan Hu (Wolfram Research) and others
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Author:  Yifan Hu              2005-04-06

#include "IpTAUCSSolverInterface.hpp"

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
# endif
#endif

namespace Ipopt
{
#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  TAUCSSolverInterface::TAUCSSolverInterface()
      :
      a_(NULL),
      multi_frontal_(false),
      negevals_(-1),
      initialized_(false)
  {
    DBG_START_METH("TAUCSSolverInterface::TAUCSSolverInterface()",dbg_verbosity);
  }

  void TAUCSSolverInterface::taucs_factor_delete_L(taucs_factorization *F)
  {
    if (F->type == TAUCS_FACTORTYPE_LLT_SUPERNODAL)
      taucs_supernodal_factor_free(F->L);
    else if (F->type == TAUCS_FACTORTYPE_LLT_CCS)
      taucs_ccs_free((taucs_ccs_matrix*) F->L);
    else if (F->type == TAUCS_FACTORTYPE_IND)
      taucs_supernodal_factor_ldlt_free(F->L);
    else
      DBG_ASSERT(false && "Factor is not deleted!");
    F->L = NULL;
  }

  void
  TAUCSSolverInterface::taucs_delete(taucs_factorization *F,
                                     taucs_ccs_matrix *A)
  {
    if (F) {
      //printf("MEM: freeing factor and colptr and rowptr\n");
      if(F->L) {
        taucs_factor_delete_L(F);
      }
      taucs_free(F->rowperm);
      taucs_free(F->colperm);
      taucs_free(F);
      F = NULL;

    }
    //printf("MEM: freeing A\n");
    if (A) {
      taucs_ccs_free(A);
      A = NULL;
    }
  }

  TAUCSSolverInterface::~TAUCSSolverInterface()
  {
    DBG_START_METH("TAUCSSolverInterface::~TAUCSSolverInterface()",
                   dbg_verbosity);

    // Tell TAUCS to release all memory
    if (initialized_) {
      taucs_delete(taucs_factor_, A_);
    }
  }

  bool TAUCSSolverInterface::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    // Tell TAUCS to release all memory if it had been used before
    if (initialized_) {
      taucs_delete(taucs_factor_, A_);
    }

    // Reset all private data
    n_ = 0;
    nz_ =0;
    initialized_ =false;

    multi_frontal_ = false;     /* we try no multifront for now */

    A_ = NULL;

    taucs_factor_ = NULL;

    return true;
  }

  ESymSolverStatus TAUCSSolverInterface::MultiSolve(bool new_matrix,
      const Index* ia,
      const Index* ja,
      Index nrhs,
      double* rhs_vals,
      bool check_NegEVals,
      Index numberOfNegEVals)
  {
    DBG_START_METH("TAUCSSolverInterface::MultiSolve",dbg_verbosity);
    DBG_ASSERT(!check_NegEVals || ProvidesInertia());
    DBG_ASSERT(initialized_);

    // check if a factorization has to be done
    if (new_matrix) {
      //delete old factor
      //printf("MEM: new matrix delete L\n");
      if (taucs_factor_->L)
        taucs_factor_delete_L(taucs_factor_);

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

  double* TAUCSSolverInterface::GetValuesArrayPtr()
  {
    DBG_ASSERT(initialized_);
    DBG_ASSERT(a_);
    return a_;
  }

  /* Initialize the local copy of the positions of the nonzero
   *  elements */
  ESymSolverStatus TAUCSSolverInterface::InitializeStructure
  (Index dim, Index nonzeros,
   const Index* ia,
   const Index* ja)
  {
    DBG_START_METH("TAUCSSolverInterface::InitializeStructure",dbg_verbosity);
    n_ = dim;
    nz_ = nonzeros;

    // create TAUCS's matrix
    //      A = (taucs_ccs_matrix*) taucs_malloc(sizeof(taucs_ccs_matrix));
    //printf("MEM: creating ccs matrix A\n");
    A_ = taucs_ccs_create(n_, n_, nz_,
                          TAUCS_SYMMETRIC | TAUCS_DOUBLE | TAUCS_LOWER);
    a_ = A_->values.d;
    for (int i = 0; i <= n_; i++)
      A_->colptr[i] = (Index) ia[i];
    for (int i = 0; i < nz_; i++)
      A_->rowind[i] = (Index) ja[i];

    // Do the symbolic facotrization
    ESymSolverStatus retval = SymbolicFactorization(ia, ja);
    if (retval != SYMSOLVER_SUCCESS) {
      return retval;
    }

    initialized_ = true;

    return retval;
  }

  ESymSolverStatus
  TAUCSSolverInterface::SymbolicFactorization(const Index* ia,
      const Index* ja)
  {
    double opt_maxdepth  = 0.0; /* default meaning no limit */
    int* rowperm = NULL;
    int* colperm = NULL;
    taucs_ccs_matrix* PAPT = NULL;

    ESymSolverStatus retcode = SYMSOLVER_SUCCESS;
    DBG_START_METH("TAUCSSolverInterface::SymbolicFactorization",
                   dbg_verbosity);

    //generate a TAUCS matrix

    //taucs factorization structure
    //printf("MEM: allocing taucs_factor\n");
    taucs_factor_ = (taucs_factorization*) taucs_malloc(sizeof(taucs_factorization));
    if (!taucs_factor_) {
      char msg[] = "MEM: taucs_factor: memory allocation\n";
      taucs_printf(msg);
      retcode = SYMSOLVER_FATAL_ERROR;
      goto release_and_return;
    }
    taucs_factor_->n       = A_->n;
    taucs_factor_->type    = TAUCS_FACTORTYPE_NONE;
    taucs_factor_->flags   = A_->flags; /* remember data type */

    // Call TAUCS to do the analysis phase
    //printf("MEM: rowperm, colperm\n");
    {
      char msg[] = "metis";
      taucs_ccs_order(A_,&rowperm,&colperm, msg);
    }
    if (!rowperm) {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "Error during ordering in TAUCS analysis phase - rowperm is NULL.\n");
      retcode = SYMSOLVER_FATAL_ERROR;
      goto release_and_return;
    }
    taucs_factor_->rowperm = rowperm;
    taucs_factor_->colperm = colperm;

    /* permute the matrix */
    //printf("MEM: making PAPT\n");
    PAPT = taucs_ccs_permute_symmetrically(A_,rowperm,colperm);
    if (!PAPT) {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "Error during permutation in TAUCS analysis phase.  PAPT is NULL.\n");
      retcode = SYMSOLVER_FATAL_ERROR;
      goto release_and_return;
    }

    // symbolic factorization
    if (multi_frontal_) {
      taucs_factor_->L = taucs_ccs_factor_ldlt_symbolic_maxdepth(PAPT,(int) opt_maxdepth);
      if (!(taucs_factor_->L)) {
        Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                       "Error during permutation in TAUCS analysis phase.  taucs_factor_->L is NULL");
        retcode = SYMSOLVER_FATAL_ERROR;
        goto release_and_return;
      }
    }
    else {
      /* for this case TAUCS does not provide separate symbolic and numerical
         factorization routines */
    }

    taucs_factor_->type = TAUCS_FACTORTYPE_IND;
    taucs_factor_->L = NULL;

    //printf("MEM: freeing PAPT\n");
    taucs_ccs_free(PAPT);
    return retcode;
release_and_return:
                      taucs_free(rowperm);
    taucs_free(colperm);
    taucs_ccs_free(PAPT);
    return retcode;
  }

  ESymSolverStatus
  TAUCSSolverInterface::Factorization(const Index* ia,
                                      const Index* ja,
                                      bool check_NegEVals,
                                      Index numberOfNegEVals)
  {
    taucs_ccs_matrix *PAPT = NULL;

    DBG_START_METH("TAUCSSolverInterface::Factorization",dbg_verbosity);

    PAPT = taucs_ccs_permute_symmetrically(A_,taucs_factor_->rowperm,
                                           taucs_factor_->colperm);

    /* numerical factorization */
    int rc = 0;
    if (multi_frontal_) {
      rc = taucs_ccs_factor_ldlt_numeric(PAPT, taucs_factor_->L);
    }
    else {
      DBG_ASSERT(!taucs_factor_->L);
      double opt_maxdepth  = 0.0; /* default meaning no limit */
      taucs_factor_->L = taucs_ccs_factor_ldlt_ll_maxdepth(PAPT,(int) opt_maxdepth);
    }
    taucs_ccs_free(PAPT);

    if (!(taucs_factor_->L) || rc) {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "Error TAUCS factorizarion phase.  rc != 0 or taucs_factor_->L is NULL\n");
      return SYMSOLVER_FATAL_ERROR;
    }


    int inertia[3];
    for (int i = 0; i < 3; i++)
      inertia[i] = 0;
    taucs_inertia_calc(taucs_factor_->L,inertia);
    negevals_ = (Index) inertia[2];


    // Check whether the number of negative eigenvalues matches the requested
    // count: what do we do with TAUCS????
    if (check_NegEVals && (numberOfNegEVals!=negevals_)) {
      return SYMSOLVER_WRONG_INERTIA;
    }

    return SYMSOLVER_SUCCESS;
  }

  ESymSolverStatus TAUCSSolverInterface::Solve(const Index* ia,
      const Index* ja,
      Index nrhs,
      double *rhs_vals)
  {
    int j;
    DBG_START_METH("TAUCSSolverInterface::Solve",dbg_verbosity);

    //printf("MEM: allocing px and pb\n");
    void *PX = (void*) taucs_malloc(element_size(A_->flags)*nrhs*(A_->n));
    void *PB = (void*) taucs_malloc(element_size(A_->flags)*nrhs*(A_->n));


    int ld = (A_->n) * element_size(A_->flags);
    for (j=0; j<nrhs; j++)
      taucs_vec_permute (A_->n, A_->flags,(char*)rhs_vals+j*ld,(char*)PB+j*ld,
                         taucs_factor_->rowperm);

    taucs_supernodal_solve_ldlt_many(taucs_factor_->L, nrhs, PX, A_->n, PB, ld);

    for (j=0; j<nrhs; j++) {
      taucs_vec_ipermute(A_->n, A_->flags, (char*)PX+j*ld, (char*) rhs_vals+j*ld,
                         taucs_factor_->rowperm);
    }

    //printf("MEM: freeing pb and px\n");
    taucs_free(PB);
    taucs_free(PX);

    return SYMSOLVER_SUCCESS;
  }

  Index TAUCSSolverInterface::NumberOfNegEVals() const
  {
    DBG_START_METH("TAUCSSolverInterface::NumberOfNegEVals",dbg_verbosity);
    DBG_ASSERT(negevals_>=0);
    return negevals_;
  }

  bool TAUCSSolverInterface::IncreaseQuality()
  {
    // At the moment, I don't see how we could tell TAUCS to do better
    return false;
  }

} // namespace Ipopt
