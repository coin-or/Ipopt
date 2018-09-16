// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Michael Hagemann               Univ of Basel 2005-10-28
//               original version (based on MA27TSolverInterface.hpp)

#ifndef __IPMA57TSOLVERINTERFACE_HPP__
#define __IPMA57TSOLVERINTERFACE_HPP__

#include "IpSparseSymLinearSolverInterface.hpp"

#ifdef FUNNY_MA57_FINT
#include <cstddef>
typedef ptrdiff_t ma57int;
#else
typedef ipfint ma57int;
#endif

namespace Ipopt
{
/** Interface to the symmetric linear solver MA57, derived from
 *  SparseSymLinearSolverInterface.
 */
class Ma57TSolverInterface: public SparseSymLinearSolverInterface
{
public:
   /** @name Constructor/Destructor */
   //@{
   /** Constructor */
   Ma57TSolverInterface();

   /** Destructor */
   virtual ~Ma57TSolverInterface();
   //@}

   bool InitializeImpl(
      const OptionsList& options,
      const std::string& prefix
      );

   /** @name Methods for requesting solution of the linear system. */
   //@{
   virtual ESymSolverStatus InitializeStructure(
      Index        dim,
      Index        nonzeros,
      const Index* airn,
      const Index* ajcn
      );

   virtual double* GetValuesArrayPtr();

   virtual ESymSolverStatus MultiSolve(
      bool         new_matrix,
      const Index* airn,
      const Index* ajcn,
      Index        nrhs,
      double*      rhs_vals,
      bool         check_NegEVals,
      Index        numberOfNegEVals
      );

   virtual Index NumberOfNegEVals() const;
   //@}

   //* @name Options of Linear solver */
   //@{
   virtual bool IncreaseQuality();

   virtual bool ProvidesInertia() const
   {
      return true;
   }

   EMatrixFormat MatrixFormat() const
   {
      return Triplet_Format;
   }
   //@}

   //@{
   static void RegisterOptions(
      SmartPtr<RegisteredOptions> roptions
      );
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
   Ma57TSolverInterface(
      const Ma57TSolverInterface&
      );

   /** Default Assignment Operator */
   void operator=(
      const Ma57TSolverInterface&
      );
   //@}

   /** @name Information about the matrix */
   //@{
   /** Number of rows and columns of the matrix */
   Index dim_;

   /** Number of nonzeros of the matrix */
   Index nonzeros_;
   //@}

   /** @name Information about most recent factorization/solve */
   //@{
   /** Number of negative eigenvalues */
   Index negevals_;
   //@}

   /** @name Initialization flags */
   //@{
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
   //@}

   /** @name Solver specific data/options */
   //@{
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
   //@}

   /** @name Data for the linear solver.
    * Storing factorization and other solver specific data structure.
    */
   //@{
   double wd_cntl_[5];
   ma57int wd_icntl_[20];

   ma57int wd_info_[40];
   double wd_rinfo_[20];

   ma57int wd_lkeep_; /* LKEEP >= 5*N + NE + max(N,NE) + 42. */
   ma57int* wd_keep_;

   ma57int* wd_iwork_; /* 5 * N. */

   double* wd_fact_;
   ma57int wd_lfact_;
   ma57int* wd_ifact_;
   ma57int wd_lifact_;

   /** factor A of matrix */
   double* a_;
   //@}

   /** @name Internal functions */
   //@{
   /** Call MA57AD and reserve memory for MA57 data.
    *
    *  Reserve memory for iw_ and ikeep_, call MA57AD to perform
    *  symbolic manipulations, and reserve all the remaining data memory
    */
   ESymSolverStatus SymbolicFactorization(
      const Index* airn,
      const Index* ajcn
      );

   /** Call MA57BD to factorize the Matrix.
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

   /** Call MA57CD to do the backsolve. */
   ESymSolverStatus Backsolve(
      Index   nrhs,
      double* rhs_vals
      );
   //@}
};

} // namespace Ipopt
#endif
