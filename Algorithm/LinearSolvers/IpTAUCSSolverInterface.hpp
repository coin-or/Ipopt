// Copyright (C) 2005, Yifan Hu (Wolfram Research) and others
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Author:  Yifan Hu              2005-04-06

#ifndef __IPTAUCSSOLVERINTERFACE_HPP__
#define __IPTAUCSSOLVERINTERFACE_HPP__

#include "IpSparseSymLinearSolverInterface.hpp"

#ifdef HAVE_CSTDLIB
# include <cstdlib>
#else
# ifdef HAVE_STDLIB_H
#  include <stdlib.h>
# else
#  error "don't have header file for stdlib"
# endif
#endif

extern "C"
{
# include <taucs.h>

  void  taucs_free_stub(void *ptr);
  void  taucs_free   (void* ptr);
# define taucs_free(x)      taucs_free_stub(x)

  void* taucs_malloc_stub(size_t size);
  void* taucs_malloc (size_t size);
# define taucs_malloc(x)    taucs_malloc_stub(x);
}

namespace Ipopt
{

  /** Interface to the linear solver TAUCS, derived from
   *  SparseSymLinearSolverInterface.  For details, see description of
   *  SparseSymLinearSolverInterface base class.
   */
  class TAUCSSolverInterface: public SparseSymLinearSolverInterface
  {
    /** Enum to define different types of factorizations that TAUCS
     *  can do. */
    enum TAUCS_Matrix_type {
      TAUCS_FACTORTYPE_NONE,
      TAUCS_FACTORTYPE_LLT_SUPERNODAL,
      TAUCS_FACTORTYPE_LLT_CCS,
      TAUCS_FACTORTYPE_IND
    };

    /** Structure to store information about the factorization */
    typedef struct
    {
      int n;
      int flags;
      TAUCS_Matrix_type type;
      int* rowperm;
      int* colperm;
      void* L;
    }
    taucs_factorization;

  public:
    /** @name Constructor/Destructor */
    //@{
    /** Constructor */
    TAUCSSolverInterface();

    /** Destructor */
    virtual ~TAUCSSolverInterface();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    bool InitializeImpl(const OptionsList& options,
                        const std::string& prefix);


    /** @name Methods for requesting solution of the linear system. */
    //@{
    /** Method for initializing internal stuctures. */
    virtual ESymSolverStatus InitializeStructure(Index dim, Index nonzeros,
        const Index *ia,
        const Index *ja);

    /** Method returing an internal array into which the nonzero
     *  elements are to be stored. */
    virtual double* GetValuesArrayPtr();

    /** Solve operation for multiple right hand sides. */
    virtual ESymSolverStatus MultiSolve(bool new_matrix,
                                        const Index* ia,
                                        const Index* ja,
                                        Index nrhs,
                                        double* rhs_vals,
                                        bool check_NegEVals,
                                        Index numberOfNegEVals);

    /** Number of negative eigenvalues detected during last
     *  factorization.
     */
    virtual Index NumberOfNegEVals() const;
    //@}

    //* @name Options of Linear solver */
    //@{
    /** Request to increase quality of solution for next solve.
     */
    virtual bool IncreaseQuality();

    /** Query whether inertia is computed by linear solver.
     *  Returns true, if linear solver provides inertia.
     */
    virtual bool ProvidesInertia() const
    {
      return true;
    }
    /** Query of requested matrix type that the linear solver
     *  understands.
     */
    EMatrixFormat MatrixFormat() const
    {
      return CSR_Format_0_Offset;
    }
    //@}

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Copy Constructor */
    TAUCSSolverInterface(const TAUCSSolverInterface&);

    /** Overloaded Equals Operator */
    void operator=(const TAUCSSolverInterface&);
    //@}

    /** @name Information about the matrix */
    //@{
    /** Number of rows and columns of the matrix */
    Index n_;

    /** Number of nonzeros of the matrix in triplet representation. */
    Index nz_;

    /** Array for storing the values of the matrix. */
    double* a_;

    /** Flag indicatining whether to use multifrontal or IC indefinite
     *  LL factorization */
    bool multi_frontal_;

    /** Structure storing information about the factorization. */
    taucs_factorization *taucs_factor_;

    /** Storing the matrix in TAUCS format */
    taucs_ccs_matrix *A_;
    //@}

    /** @name Information about most recent factorization/solve */
    //@{
    /** Number of negative eigenvalues */
    Index negevals_;

    /** Flag indicating if internal data is initialized.
     *  For initialization, this object needs to have seen a matrix */
    bool initialized_;
    //@}

    /** @name Auxilliary methods */
    //@{
    /** TAUCS to do the analysis phase. */
    ESymSolverStatus SymbolicFactorization(const Index* ia,
                                           const Index* ja);

    /** Call TAUCS to factorize the Matrix. */
    ESymSolverStatus Factorization(const Index* ia,
                                   const Index* ja,
                                   bool check_NegEVals,
                                   Index numberOfNegEVals);

    /** Call TAUCS to do the Solve. */
    ESymSolverStatus Solve(const Index* ia,
                           const Index* ja,
                           Index nrhs,
                           double *rhs_vals);

    /** return size of matrix entry type */
    int element_size(int flags)
    {
      if (flags & TAUCS_SINGLE)
        return sizeof(taucs_single);
      if (flags & TAUCS_DOUBLE)
        return sizeof(taucs_double);
      if (flags & TAUCS_SCOMPLEX)
        return sizeof(taucs_scomplex);
      if (flags & TAUCS_DCOMPLEX)
        return sizeof(taucs_dcomplex);
      if (flags & TAUCS_INT)
        return sizeof(int);
      return -1;
    }

    /** Delete TAUCS' factorization structure */
    void taucs_factor_delete_L(taucs_factorization *F);

    /** Delete TAUCS' factorization and matrix */
    void taucs_delete(taucs_factorization *F, taucs_ccs_matrix *A);
    //@}
  };

} // namespace Ipopt
#endif
