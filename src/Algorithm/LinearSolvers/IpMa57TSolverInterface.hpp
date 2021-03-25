// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Michael Hagemann               Univ of Basel 2005-10-28
//               original version (based on MA27TSolverInterface.hpp)

#ifndef __IPMA57TSOLVERINTERFACE_HPP__
#define __IPMA57TSOLVERINTERFACE_HPP__

#include "IpSparseSymLinearSolverInterface.hpp"
#include "IpLibraryLoader.hpp"
#include "IpTypes.h"

#ifdef FUNNY_MA57_FINT
#include <cstddef>
typedef ptrdiff_t ma57int;
#else
#include "IpTypes.h"
typedef ipfint ma57int;
#endif

#define IPOPT_DECL_MA57A(x) void (x)( \
   ipfint*       n,     /**< Order of matrix. */ \
   ipfint*       ne,    /**< Number of entries. */ \
   const ipfint* irn,   /**< Matrix nonzero row structure */ \
   const ipfint* jcn,   /**< Matrix nonzero column structure */ \
   ipfint*       lkeep, /**< Workspace for the pivot order of lenght 3*n */ \
   ipfint*       keep,  /**< Workspace for the pivot order of lenght 3*n */ \
   /* Automatically iflag = 0; ikeep pivot order iflag = 1 */ \
   ipfint*       iwork, /**< Integer work space. */ \
   ipfint*       icntl, /**< Integer Control parameter of length 30 */ \
   ipfint*       info,  /**< Statistical Information; Integer array of length 20 */ \
   ipnumber*     rinfo  /**< Double Control parameter of length 5 */ \
)

#define IPOPT_DECL_MA57B(x) void (x)( \
   ipfint*   n,      /**< Order of matrix. */ \
   ipfint*   ne,     /**< Number of entries. */ \
   ipnumber* a,      /**< Numerical values. */ \
   ipnumber* fact,   /**< Entries of factors. */ \
   ipfint*   lfact,  /**< Length of array `fact'. */ \
   ipfint*   ifact,  /**< Indexing info for factors. */ \
   ipfint*   lifact, /**< Length of array `ifact'. */ \
   ipfint*   lkeep,  /**< Length of array `keep'. */ \
   ipfint*   keep,   /**< Integer array. */ \
   ipfint*   iwork,  /**< Workspace of length `n'. */ \
   ipfint*   icntl,  /**< Integer Control parameter of length 20. */ \
   ipnumber* cntl,   /**< Double Control parameter of length 5. */ \
   ipfint*   info,   /**< Statistical Information; Integer array of length 40. */ \
   ipnumber* rinfo   /**< Statistical Information; Real array of length 20. */ \
)

/* Solution job:  Solve for...
 * - JOB <= 1:  A
 * - JOB == 2:  PLP^t
 * - JOB == 3:  PDP^t
 * - JOB >= 4:  PL^t P^t
 */
#define IPOPT_DECL_MA57C(x) void (x)( \
   ipfint*   job,    /**< Solution job. */ \
   ipfint*   n,      /**< Order of matrix. */ \
   ipnumber* fact,   /**< Entries of factors. */ \
   ipfint*   lfact,  /**< Length of array `fact'. */ \
   ipfint*   ifact,  /**< Indexing info for factors. */ \
   ipfint*   lifact, /**< Length of array `ifact'. */ \
   ipfint*   nrhs,   /**< Number of right hand sides. */ \
   ipnumber* rhs,    /**< Numerical Values. */ \
   ipfint*   lrhs,   /**< Leading dimensions of `rhs'. */ \
   ipnumber* work,   /**< Real workspace. */ \
   ipfint*   lwork,  /**< Length of `work', >= N*NRHS. */ \
   ipfint*   iwork,  /**< Integer array of length `n'. */ \
   ipfint*   icntl,  /**< Integer Control parameter array of length 20. */ \
   ipfint*   info    /**< Statistical Information; Integer array of length 40. */ \
)

#define IPOPT_DECL_MA57E(x) void (x)( \
   ipfint*   n, \
   ipfint*   ic,    /**< 0: copy real array.  >=1:  copy integer array. */ \
   ipfint*   keep,   \
   ipnumber* fact,   \
   ipfint*   lfact,  \
   ipnumber* newfac, \
   ipfint*   lnew,   \
   ipfint*   ifact,  \
   ipfint*   lifact, \
   ipfint*   newifc, \
   ipfint*   linew,  \
   ipfint*   info    \
)

#define IPOPT_DECL_MA57I(x) void (x)( \
   ipnumber* cntl, \
   ipfint*   icntl \
)

namespace Ipopt
{
/** Interface to the symmetric linear solver MA57, derived from
 *  SparseSymLinearSolverInterface.
 */
class Ma57TSolverInterface: public SparseSymLinearSolverInterface
{
public:
   /** @name Constructor/Destructor */
   ///@{
   /** Constructor */
   Ma57TSolverInterface(
      SmartPtr<LibraryLoader> hslloader_
   );

   /** Destructor */
   virtual ~Ma57TSolverInterface();
   ///@}

   bool InitializeImpl(
      const OptionsList& options,
      const std::string& prefix
   );

   /** @name Methods for requesting solution of the linear system. */
   ///@{
   virtual ESymSolverStatus InitializeStructure(
      Index        dim,
      Index        nonzeros,
      const Index* airn,
      const Index* ajcn
   );

   virtual Number* GetValuesArrayPtr();

   virtual ESymSolverStatus MultiSolve(
      bool         new_matrix,
      const Index* airn,
      const Index* ajcn,
      Index        nrhs,
      Number*      rhs_vals,
      bool         check_NegEVals,
      Index        numberOfNegEVals
   );

   virtual Index NumberOfNegEVals() const;
   ///@}

   //* @name Options of Linear solver */
   ///@{
   virtual bool IncreaseQuality();

   virtual bool ProvidesInertia() const
   {
      return true;
   }

   EMatrixFormat MatrixFormat() const
   {
      return Triplet_Format;
   }
   ///@}

   ///@{
   static void RegisterOptions(
      SmartPtr<RegisteredOptions> roptions
   );
   ///@}

   /// set MA57 functions to use for every instantiation of this class
   static void SetFunctions(
      IPOPT_DECL_MA57A(*ma57a),
      IPOPT_DECL_MA57B(*ma57b),
      IPOPT_DECL_MA57C(*ma57c),
      IPOPT_DECL_MA57E(*ma57e),
      IPOPT_DECL_MA57I(*ma57i)
      );

private:
   /**@name Default Compiler Generated Methods
    * (Hidden to avoid implicit creation/calling).
    * These methods are not implemented and
    * we do not want the compiler to implement
    * them for us, so we declare them private
    * and do not define them. This ensures that
    * they will not be implicitly created/called. */
   ///@{
   /** Copy Constructor */
   Ma57TSolverInterface(
      const Ma57TSolverInterface&
   );

   /** Default Assignment Operator */
   void operator=(
      const Ma57TSolverInterface&
   );
   ///@}

   /**@name MA57 function pointers
    * @{
    */
   SmartPtr<LibraryLoader> hslloader;

   /// symbolic factorization
   IPOPT_DECL_MA57A(*ma57a);
   /// numerical factorization
   IPOPT_DECL_MA57B(*ma57b);
   /// solution
   IPOPT_DECL_MA57C(*ma57c);
   /// copy arrays
   IPOPT_DECL_MA57E(*ma57e);
   /// initialize solver
   IPOPT_DECL_MA57I(*ma57i);
   ///@}

   /** @name Information about the matrix */
   ///@{
   /** Number of rows and columns of the matrix */
   Index dim_;

   /** Number of nonzeros of the matrix */
   Index nonzeros_;
   ///@}

   /** @name Information about most recent factorization/solve */
   ///@{
   /** Number of negative eigenvalues */
   Index negevals_;
   ///@}

   /** @name Initialization flags */
   ///@{
   /** Flag indicating if internal data is initialized.
    *
    *  For initialization, this object needs to have seen a matrix
    */
   bool initialized_;
   /** Flag indicating if the matrix has to be refactorized because
    *  the pivot tolerance has been changed.
    */
   bool pivtol_changed_;
   /** Flag that is true if we just requested the values of the
    *  matrix again (SYMSOLVER_CALL_AGAIN) and have to factorize
    *  again.
    */
   bool refactorize_;
   ///@}

   /** @name Solver specific data/options */
   ///@{
   /** Pivot tolerance */
   Number pivtol_;
   /** Maximal pivot tolerance */
   Number pivtolmax_;
   /** Factor for estimating initial size of work arrays */
   Number ma57_pre_alloc_;
   /** Flag indicating whether the TNLP with identical structure has
    *  already been solved before.
    */
   bool warm_start_same_structure_;
   ///@}

   /** @name Data for the linear solver.
    * Storing factorization and other solver specific data structure.
    */
   ///@{
   Number wd_cntl_[5];
   ma57int wd_icntl_[20];

   ma57int wd_info_[40];
   Number wd_rinfo_[20];

   ma57int wd_lkeep_; /* LKEEP >= 5*N + NE + max(N,NE) + 42. */
   ma57int* wd_keep_;

   ma57int* wd_iwork_; /* 5 * N. */

   Number* wd_fact_;
   ma57int wd_lfact_;
   ma57int* wd_ifact_;
   ma57int wd_lifact_;

   /** factor A of matrix */
   Number* a_;
   ///@}

   /** @name Internal functions */
   ///@{
   /** Call MA57AX and reserve memory for MA57 data.
    *
    *  Reserve memory for iw_ and ikeep_, call MA57AD to perform
    *  symbolic manipulations, and reserve all the remaining data memory
    */
   ESymSolverStatus SymbolicFactorization(
      const Index* airn,
      const Index* ajcn
   );

   /** Call MA57BX to factorize the Matrix.
    *
    *  It is assumed that the first nonzeros_ element of a_ contain the values
    *  of the matrix to be factorized.
    */
   ESymSolverStatus Factorization(
      const Index* airn,
      const Index* ajcn,
      bool         check_NegEVals,
      Index        numberOfNegEVals
   );

   /** Call MA57CX to do the backsolve. */
   ESymSolverStatus Backsolve(
      Index   nrhs,
      Number* rhs_vals
   );
   ///@}
};

} // namespace Ipopt
#endif
