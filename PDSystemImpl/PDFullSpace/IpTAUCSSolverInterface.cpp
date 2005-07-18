// Copyright (C) 2005, Yifan Hu (Wolfram Research)
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Author:  Yifan Hu              2005-04-06

#include "IpTAUCSSolverInterface.hpp"

#ifdef OLD_C_HEADERS
# include <math.h>
#else
# include <cmath>
#endif

namespace Ipopt
{

  DBG_SET_VERBOSITY(0);

  TAUCSSolverInterface::TAUCSSolverInterface()
      :
      m(0),
      n(0),
      nz(0),
      initialized(false),
      negevals(-1),
      a(NULL),
      multi_frontal(false)
  {
    DBG_START_METH("TAUCSSolverInterface::TAUCSSolverInterface()",dbg_verbosity);

  }

  static void taucs_factor_delete_L(taucs_factorization *F)
  {
    if (F->type == TAUCS_FACTORTYPE_LLT_SUPERNODAL)
      taucs_supernodal_factor_free(F->L);
    if (F->type == TAUCS_FACTORTYPE_LLT_CCS)
      taucs_ccs_free((taucs_ccs_matrix*) F->L);
  }

  static void taucs_delete(taucs_factorization *F, taucs_ccs_matrix *A)
  {
    if (F) {
      //printf("MEM: freeing factor and colptr and rowptr\n");
      taucs_factor_delete_L(F);
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
    if (initialized) {
      taucs_delete(taucs_factor, A);
      //         DBG_ASSERT(ERROR==0);
    }
  }

  bool TAUCSSolverInterface::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    Number value = 0.0;

    // Tell TAUCS to release all memory if it had been used before
    if (initialized) {
      taucs_delete(taucs_factor, A);
      //   DBG_ASSERT(ERROR==0);
    }

    // Reset all private data
    m = n = 0;
    nz =0;
    initialized =false;

    multi_frontal = false;     /* we try no multifront for now */

    A = NULL;

    taucs_factor = NULL;

    // Call TAUCS's initialization routine

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
    DBG_ASSERT(initialized);

    // check if a factorization has to be done
    if (new_matrix) {
      //delete old factor
      //printf("MEM: new matrix delete L\n");
      taucs_factor_delete_L(taucs_factor);

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
    DBG_ASSERT(initialized);
    DBG_ASSERT(a);
    return a;
  }

  /** Initialize the local copy of the positions of the nonzero
      elements */
  ESymSolverStatus TAUCSSolverInterface::InitializeStructure
  (Index dim, Index nonzeros,
   const Index* ia,
   const Index* ja)
  {
    DBG_START_METH("TAUCSSolverInterface::InitializeStructure",dbg_verbosity);
    m = n = dim;
    nz = nonzeros;

    // create TAUCS's matrix
    //      A = (taucs_ccs_matrix*) taucs_malloc(sizeof(taucs_ccs_matrix));
    //printf("MEM: creating ccs matrix A\n");
    A = taucs_ccs_create(m, n, nz, TAUCS_SYMMETRIC | TAUCS_DOUBLE | TAUCS_LOWER);
    a = A->values.d;
    for (int i = 0; i <= m; i++)
      A->colptr[i] = (Index) ia[i];
    for (int i = 0; i < nz; i++)
      A->rowind[i] = (Index) ja[i];

    // Do the symbolic facotrization
    ESymSolverStatus retval = SymbolicFactorization(ia, ja);
    if (retval != SYMSOLVER_SUCCESS) {
      return retval;
    }

    initialized = true;

    return retval;
  }

  ESymSolverStatus
  TAUCSSolverInterface::SymbolicFactorization(const Index* ia,
      const Index* ja)
  {
    double opt_maxdepth  = 0.0; /* default meaning no limit */
    int* rowperm = NULL, *colperm = NULL;
    taucs_ccs_matrix *PAPT = NULL;

    ESymSolverStatus retcode = SYMSOLVER_SUCCESS;
    DBG_START_METH("TAUCSSolverInterface::SymbolicFactorization",
                   dbg_verbosity);

    //generate a TAUCS matrix

    //taucs factorization structure
    //printf("MEM: allocing taucs_factor\n");
    taucs_factor = (taucs_factorization*) taucs_malloc(sizeof(taucs_factorization));
    if (!taucs_factor) {
      taucs_printf("MEM: taucs_factor: memory allocation\n");
      retcode = SYMSOLVER_FATAL_ERROR;
      goto release_and_return;
    }
    taucs_factor->n       = A->n;
    taucs_factor->type    = TAUCS_FACTORTYPE_NONE;
    taucs_factor->flags   = A->flags; /* remember data type */

    // Call TAUCS to do the analysis phase
    //printf("MEM: rowperm, colperm\n");
    taucs_ccs_order(A,&rowperm,&colperm, "metis");
    if (!rowperm) {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "Error during ordering in TAUCS analysis phase.  ERROR = %d.\n",
                     SYMSOLVER_FATAL_ERROR);
      retcode = SYMSOLVER_FATAL_ERROR;
      goto release_and_return;
    }
    taucs_factor->rowperm = rowperm;
    taucs_factor->colperm = colperm;

    /* permute the matrix */
    //printf("MEM: making PAPT\n");
    PAPT = taucs_ccs_permute_symmetrically(A,rowperm,colperm);
    if (!PAPT) {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "Error during permutation in TAUCS analysis phase.  ERROR = %d.\n",
                     SYMSOLVER_FATAL_ERROR);
      retcode = SYMSOLVER_FATAL_ERROR;
      goto release_and_return;
    }

    // symbolic factorization
    if (multi_frontal) {
      taucs_factor->L = taucs_ccs_factor_ldlt_symbolic_maxdepth(PAPT,(int) opt_maxdepth);
      if (!(taucs_factor->L)) {
        Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                       "Error during permutation in TAUCS analysis phase.  ERROR = %d.\n",
                       SYMSOLVER_FATAL_ERROR);
        retcode = SYMSOLVER_FATAL_ERROR;
        goto release_and_return;
      }
    }
    else {
      /* for this case TAUCS does not provide separate symbolic and numerical
         factorization routines */
    }

    taucs_factor->type = TAUCS_FACTORTYPE_IND;

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

    PAPT = taucs_ccs_permute_symmetrically(A,taucs_factor->rowperm,
                                           taucs_factor->colperm);

    /* numerical factorization */
    int rc = 0;
    if (multi_frontal) {
      rc = taucs_ccs_factor_ldlt_numeric(PAPT, taucs_factor->L);
    }
    else {

      double opt_maxdepth  = 0.0; /* default meaning no limit */
      taucs_factor->L = taucs_ccs_factor_ldlt_ll_maxdepth(PAPT,(int) opt_maxdepth);

    }
    taucs_ccs_free(PAPT);

    if (!(taucs_factor->L) || rc) {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "Error TAUCS factorizarion phase.  ERROR = %d.\n",
                     SYMSOLVER_FATAL_ERROR);
      return SYMSOLVER_FATAL_ERROR;
    }


    int inertia[3];
    for (int i = 0; i < 3; i++)
      inertia[i] = 0;
    taucs_inertia_calc(taucs_factor->L,inertia);
    negevals = (Index) inertia[2];


    // Check whether the number of negative eigenvalues matches the requested
    // count: what do we do with TAUCS????
    if (check_NegEVals && (numberOfNegEVals!=negevals)) {
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
    void *PX = (void*) taucs_malloc(element_size(A->flags)*nrhs*(A->n));
    void *PB = (void*) taucs_malloc(element_size(A->flags)*nrhs*(A->n));


    int ld = (A->n) * element_size(A->flags);
    for (j=0; j<nrhs; j++)
      taucs_vec_permute (A->n, A->flags,(char*)rhs_vals+j*ld,(char*)PB+j*ld,
                         taucs_factor->rowperm);

    taucs_supernodal_solve_ldlt_many(taucs_factor->L, nrhs, PX, A->n, PB, ld);

    for (j=0; j<nrhs; j++) {
      taucs_vec_ipermute(A->n, A->flags, (char*)PX+j*ld, (char*) rhs_vals+j*ld,
                         taucs_factor->rowperm);
    }

    //printf("MEM: freeing pb and px\n");
    taucs_free(PB);
    taucs_free(PX);

    return SYMSOLVER_SUCCESS;
  }

  Index TAUCSSolverInterface::NumberOfNegEVals() const
  {
    DBG_START_METH("TAUCSSolverInterface::NumberOfNegEVals",dbg_verbosity);
    DBG_ASSERT(negevals>=0);
    return negevals;
  }

  bool TAUCSSolverInterface::IncreaseQuality()
  {
    // At the moment, I don't see how we could tell TAUCS to do better
    return false;
  }

} // namespace Ipopt
