// Copyright (C) 2005, Yifan Hu (Wolfram Research)
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Author:  Yifan Hu              2005-04-06

#ifndef __IPTAUCSSOLVERINTERFACE_HPP__
#define __IPTAUCSSOLVERINTERFACE_HPP__

#include "IpSparseSymLinearSolverInterface.hpp"

#include <stdlib.h>

extern "C"
{
  // we can not just wrap a extern C around
  // #include "taucs.h"
  // becuase we will get tons of messages and errors,
  // because taucs.h includes a lot of C libraried.
  typedef double    taucs_double;
  typedef float     taucs_single;
  //we do not use complex anyway, so do worry about
  // whether TAUCS_C99_COMPLEX is defined or not
#ifdef TAUCS_C99_COMPLEX
#include <complex.h>

  typedef _Complex double taucs_dcomplex;
  typedef _Complex float  taucs_scomplex;
#else /* C99 */

  typedef struct {
    double r,i;
  }
  taucs_dcomplex;
  typedef struct {
    float  r,i;
  }
  taucs_scomplex;
#endif
#define TAUCS_INT       1024
#define TAUCS_DOUBLE    2048
#define TAUCS_SINGLE    4096
#define TAUCS_DCOMPLEX  8192
#define TAUCS_SCOMPLEX 16384
#define TAUCS_LOWER      1
#define TAUCS_UPPER      2
#define TAUCS_TRIANGULAR 4
#define TAUCS_SYMMETRIC  8
#define TAUCS_HERMITIAN  16
#define TAUCS_PATTERN    32
#define TAUCS_FACTORTYPE_NONE           0

  typedef struct {
    int   n;
    int   flags;
    int   type;
    int*  rowperm;
    int*  colperm;
    void* L;
  }
  taucs_factorization;
  typedef struct {
    int     n;    /* columns                      */
    int     m;    /* rows; don't use if symmetric   */
    int     flags;
    int*    colptr; /* pointers to where columns begin in rowind and values. */
    /* 0-based. Length is (n+1). */
    int*    rowind; /* row indices */

    union {
      void*           v;
      taucs_double*   d;
      taucs_single*   s;
      taucs_dcomplex* z;
      taucs_scomplex* c;
    } values;

  }
  taucs_ccs_matrix;


  void  taucs_free_stub(void *ptr);
  void  taucs_free   (void* ptr);
#define taucs_free(x)      taucs_free_stub(x)
  /*
   extern void iFree(void *bp);
  #define taucs_free(x)  iFree(x)
  */


  void* taucs_malloc_stub(size_t size);
  void* taucs_malloc (size_t size);
#define taucs_malloc(x)    taucs_malloc_stub(x);
  /*
   extern void *MallocAlign16( size_t size);
  #define taucs_malloc(x)    MallocAlign16(x)
  */

  int taucs_printf(char *fmt, ...);
  void taucs_ccs_free(taucs_ccs_matrix* matrix);
  void taucs_ccs_order(taucs_ccs_matrix* matrix,
                       int** perm, int** invperm,
                       char* which);
  void taucs_vec_permute(int n, int flags, void* v, void* pv, int p[]);
  void    taucs_inertia_calc(void* vL, int* inertia);
  void taucs_supernodal_factor_free                (void* L);
  void taucs_supernodal_factor_free_numeric        (void* L);
  int   taucs_ccs_factor_ldlt_numeric(taucs_ccs_matrix* A,void* L);
  taucs_ccs_matrix* taucs_ccs_permute_symmetrically (taucs_ccs_matrix* A,
      int* perm, int* invperm);
  int taucs_supernodal_solve_ldlt_many(void *L,int n,void* X, int ld_X,void* B,
                                       int ld_B);
  void* taucs_ccs_factor_ldlt_ll_maxdepth(taucs_ccs_matrix* A,int max_depth);
  void* taucs_ccs_factor_ldlt_symbolic_maxdepth(taucs_ccs_matrix* A,int max_depth);
  void taucs_vec_ipermute(int n, int flags, void* v, void* pv, int p[]);
  taucs_ccs_matrix* taucs_ccs_create(int m, int n, int nnz, int flags);

}

static int element_size(int flags)
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

/*********************************************************/
/* Generic Factor routines                               */
/* (Experimental, unstable interface)                    */
/*********************************************************/

#define TAUCS_FACTORTYPE_NONE           0
#define TAUCS_FACTORTYPE_LLT_SUPERNODAL 1
#define TAUCS_FACTORTYPE_LLT_CCS        2
#define TAUCS_FACTORTYPE_LDLT_CCS       3
#define TAUCS_FACTORTYPE_LLT_OOC        4
#define TAUCS_FACTORTYPE_LU_OOC         5
#define TAUCS_FACTORTYPE_IND            6
#define TAUCS_FACTORTYPE_IND_OOC        7
#define TAUCS_FACTORTYPE_LU             8


namespace Ipopt
{

  /** Interface to the linear solver TAUCS, derived from
   *  SparseSymLinearSolverInterface.  For details, see description of
   *  SparseSymLinearSolverInterface base class.
   */
  class TAUCSSolverInterface: public SparseSymLinearSolverInterface
  {
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
    Index m, n;

    /** Number of nonzeros of the matrix in triplet representation. */
    Index nz;

    /** Array for storing the values of the matrix. */
    double* a;

    /** @name Information about most recent factorization/solve */
    /** Number of negative eigenvalues */
    Index negevals;

    /** @name Initialization flags */
    /** Flag indicating if internal data is initialized.
     *  For initialization, this object needs to have seen a matrix */
    bool initialized;

    /* whether to use multifrontal or IC indefinite LL factorization */
    bool multi_frontal;

    /* name Solver specific information */
    taucs_factorization *taucs_factor;
    taucs_ccs_matrix *A;

    // Call TAUCS to do the analysis phase.
    ESymSolverStatus SymbolicFactorization(const Index* ia,
                                           const Index* ja);

    // Call TAUCS to factorize the Matrix.
    ESymSolverStatus Factorization(const Index* ia,
                                   const Index* ja,
                                   bool check_NegEVals,
                                   Index numberOfNegEVals);

    // Call TAUCS to do the Solve.
    ESymSolverStatus Solve(const Index* ia,
                           const Index* ja,
                           Index nrhs,
                           double *rhs_vals);
    //@}
  };

} // namespace Ipopt
#endif
