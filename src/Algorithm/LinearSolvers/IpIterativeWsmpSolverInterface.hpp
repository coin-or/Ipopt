// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Andreas Waechter              IBM    2009-09-18
//               based on IpWsmpSolverInterface.hpp (rev 1483)

#ifndef __IPITERATIVEWSMPSOLVERINTERFACE_HPP__
#define __IPITERATIVEWSMPSOLVERINTERFACE_HPP__

#include "IpSparseSymLinearSolverInterface.hpp"
#include "IpTypes.h"

namespace Ipopt
{

/** Interface to the linear solver WISMP, derived from
 *  SparseSymLinearSolverInterface.
 */
class IterativeWsmpSolverInterface: public SparseSymLinearSolverInterface
{
public:
   /** @name Constructor/Destructor */
   ///@{
   /** Constructor */
   IterativeWsmpSolverInterface();

   /** Destructor */
   virtual ~IterativeWsmpSolverInterface();
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
      return false;
   }

   EMatrixFormat MatrixFormat() const
   {
      return CSR_Format_1_Offset;
   }
   ///@}

   static void RegisterOptions(
      SmartPtr<RegisteredOptions> roptions
   );

private:
   /**@name Default Compiler Generated Methods
    * (Hidden to avoid implicit creation/calling).
    * These methods are not implemented and
    * we do not want the compiler to implement
    * them for us, so we declare them private
    * and do not define them. This ensures that
    * they will not be implicitly created/called.
    */
   ///@{
   /** Copy Constructor */
   IterativeWsmpSolverInterface(
      const IterativeWsmpSolverInterface&
   );

   /** Default Assignment Operator */
   void operator=(
      const IterativeWsmpSolverInterface&
   );
   ///@}

   /** @name Information about the matrix */
   ///@{
   /** Number of rows and columns of the matrix */
   Index dim_;

   /** Array for storing the values of the matrix. */
   double* a_;
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
   /** iteration number in which matrices are to be written out */
   Index wsmp_write_matrix_iteration_;
   Number wsmp_inexact_droptol_;
   Number wsmp_inexact_fillin_limit_;
   ///@}

   /** Counter for matrix file numbers */
   Index matrix_file_number_;

   /** @name Information about most recent factorization/solve */
   ///@{
#if 0
   /** Number of negative eigenvalues */
   Index negevals_;
#endif
   ///@}

   /** @name Initialization flags */
   ///@{
   /** Flag indicating if internal data is initialized.
    *  For initialization, this object needs to have seen a matrix
    */
   bool initialized_;
   /** Flag indicating if the matrix has to be refactorized because
    *  the pivot tolerance has been changed.
    */
   bool pivtol_changed_;
   /** Flag indicating whether symbolic factorization and order has
    *  already been performed.
    */
   bool have_symbolic_factorization_;
   ///@}

   /** @name Solver specific information */
   ///@{
   /** Integer parameter array for WISMP. */
   ipindex* IPARM_;
   /** Double precision parameter array for WISMP. */
   double* DPARM_;
   ///@}

   /** @name Internal functions */
   ///@{
   /** Call Wsmp to do the analysis phase.
    */
   ESymSolverStatus SymbolicFactorization(
      const Index* ia,
      const Index* ja
   );

   /** Call Wsmp to really do the analysis phase. */
   ESymSolverStatus InternalSymFact(
      const Index* ia,
      const Index* ja
   );

   /** Call Wsmp to factorize the Matrix. */
   ESymSolverStatus Factorization(
      const Index* ia,
      const Index* ja,
      bool         check_NegEVals,
      Index        numberOfNegEVals
   );

   /** Call Wsmpx to do the Solve. */
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
