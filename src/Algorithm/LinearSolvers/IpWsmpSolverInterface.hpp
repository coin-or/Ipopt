// Copyright (C) 2005, 2010 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-03-17

#ifndef __IPWSMPSOLVERINTERFACE_HPP__
#define __IPWSMPSOLVERINTERFACE_HPP__

#include "IpSparseSymLinearSolverInterface.hpp"
#include "IpPardisoSolverInterface.hpp"  // for IPOPT_DECL_SMAT_REORDERING_PARDISO_WSMP
#include "IpTypes.h"

namespace Ipopt
{

/** Interface to the linear solver Wsmp, derived from
 *  SparseSymLinearSolverInterface.
 */
class WsmpSolverInterface: public SparseSymLinearSolverInterface
{
public:
   /** @name Constructor/Destructor */
   ///@{
   /** Constructor */
   WsmpSolverInterface(
#ifdef PARDISO_MATCHING_PREPROCESS
      SmartPtr<LibraryLoader> pardisoloader_  ///< @since 3.14.0
#endif
   );

   /** Destructor */
   virtual ~WsmpSolverInterface();
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
      const Index* ia,
      const Index* ja
   );

   virtual double* GetValuesArrayPtr();

   virtual ESymSolverStatus MultiSolve(
      bool         new_matrix,
      const Index* ia,
      const Index* ja,
      Index        nrhs,
      double*      rhs_vals,
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
      return CSR_Format_1_Offset;
   }
   ///@}

   ///@{
   static void RegisterOptions(
      SmartPtr<RegisteredOptions> roptions
   );

   /// give WSMP version
   static void GetVersion(
      int& V,
      int& R,
      int& M
   );
   ///@}

   virtual bool ProvidesDegeneracyDetection() const;

   virtual ESymSolverStatus DetermineDependentRows(
      const Index*      ia,
      const Index*      ja,
      std::list<Index>& c_deps
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
   WsmpSolverInterface(
      const WsmpSolverInterface&
   );

   /** Default Assignment Operator */
   void operator=(
      const WsmpSolverInterface&
   );
   ///@}

   /** @name Information about the matrix */
   ///@{
   /** Number of rows and columns of the matrix */
   Index dim_;

   /** Number of nonzeros of the matrix in triplet representation. */
   Index nonzeros_;

   /** Array for storing the values of the matrix. */
   double* a_;

#ifdef PARDISO_MATCHING_PREPROCESS
   /**  @name Arrays for storing the values of a second matrix that has been already reordered. */
   ///@{
   Index* ia2;
   Index* ja2;
   double* a2_;
   Index* perm2;
   double* scale2;
   ///@}
#endif

   ///@}

   /** @name Solver specific options */
   ///@{
   /** Option that controls the matching strategy. */
   Index wsmp_num_threads_;
   /** Pivot tolerance */
   Number wsmp_pivtol_;
   /** Maximal pivot tolerance */
   Number wsmp_pivtolmax_;
   /** Indicating which of WSMP's scaling methods should be used. */
   Index wsmp_scaling_;
   /** WSMP's singularity threshold.
    *
    *  The smaller this value the less likely a matrix is declared singular.
    */
   Number wsmp_singularity_threshold_;
   /** iteration number in which matrices are to be written out */
   Index wsmp_write_matrix_iteration_;
   /** Flag indicating if the inertia is always assumed to be correct. */
   bool skip_inertia_check_;
   /** Flag indicating whether the positive definite version of WSMP should be used */
   bool wsmp_no_pivoting_;
   ///@}

   /** Counter for matrix file numbers */
   Index matrix_file_number_;

   /** @name Information about most recent factorization/solve */
   ///@{
   /** Number of negative eigenvalues */
   Index negevals_;
   ///@}

   /** @name Initialization flags */
   ///@{
   /** Flag indicating if internal data is initialized.
    *
    *  For initialization, this object needs to have seen a matrix.
    */
   bool initialized_;
   /** Flag indicating if we already printed how many threads are used by WSMP. */
   bool printed_num_threads_;
   /** Flag indicating if the matrix has to be refactorized because
    *  the pivot tolerance has been changed, or the computation of
    *  the ordering has been triggered with DPARNM[14].
    */
   bool pivtol_changed_;
   /** Flag indicating whether symbolic factorization and order has
    *  already been performed.
    */
   bool have_symbolic_factorization_;
   /** Counter indicating how many factorizations have been done sine
    *  the last recomputation of the ordering.
    */
   Index factorizations_since_recomputed_ordering_;
   ///@}

   /** @name Solver specific information */
   ///@{
   /** Integer parameter array for WSSMP. */
   Index* IPARM_;
   /** Double precision parameter array for WSSMP. */
   double* DPARM_;
   /** WSSMP's permutation vector */
   Index* PERM_;
   /** WSSMP's inverse permutation vector */
   Index* INVP_;
   /** WSSMP's internal MRP array */
   Index* MRP_;
   ///@}

   /**@name PARDISO function pointer
    * @{
    */
#ifdef PARDISO_MATCHING_PREPROCESS
   SmartPtr<LibraryLoader> pardisoloader;
   IPOPT_DECL_SMAT_REORDERING_PARDISO_WSMP(*smat_reordering_pardiso_wsmp);
#endif
   /**@} */

   /** @name Internal functions */
   ///@{
   /** Call Wsmp to do the analysis phase. */
   ESymSolverStatus SymbolicFactorization(
      const Index* ia,
      const Index* ja
   );

   /** Call Wsmp to really do the analysis phase. */
   ESymSolverStatus InternalSymFact(
      const Index* ia,
      const Index* ja,
      Index        numberOfNegEVals
   );

   /** Call Wsmp to factorize the Matrix. */
   ESymSolverStatus Factorization(
      const Index* ia,
      const Index* ja,
      bool         check_NegEVals,
      Index        numberOfNegEVals
   );

   /** Call Wsmp to do the Solve. */
   ESymSolverStatus Solve(
      const Index* ia,
      const Index* ja,
      Index        nrhs,
      double*      rhs_vals
   );
   ///@}
};

} // namespace Ipopt
#endif
