// Copyright (C) 2006, 2007 Damien Hocking, KBC Advanced Technologies
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors: Damien Hocking                 KBC    2006-03-20
//        (included his original contribution into Ipopt package on 2006-03-25)
//          Andreas Waechter               IBM    2006-03-25
//           (minor changes and corrections)
//          Scott Turnberg                 CMU    2006-05-12
//           (major revision)
//           (incorporated by AW on 2006-11-11 into Ipopt package)

#ifndef __IPMUMPSSOLVERINTERFACE_HPP__
#define __IPMUMPSSOLVERINTERFACE_HPP__

#include "IpSparseSymLinearSolverInterface.hpp"

namespace Ipopt
{

/** Interface to the linear solver Mumps, derived from
 *  SparseSymLinearSolverInterface.
 */
class MumpsSolverInterface: public SparseSymLinearSolverInterface
{
public:
   /** @name Constructor/Destructor */
   ///@{
   /** Constructor */
   MumpsSolverInterface();

   /** Destructor */
   virtual ~MumpsSolverInterface();
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

   static void RegisterOptions(
      SmartPtr<RegisteredOptions> roptions
   );

   /// give name of MUMPS with version info
   /// @since 3.14.0
   static std::string GetName();

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
   MumpsSolverInterface(
      const MumpsSolverInterface&
   );

   /** Default Assignment Operator */
   void operator=(
      const MumpsSolverInterface&
   );
   ///@}

   /** @name Information about the matrix */
   ///@{
   /** Primary MUMP data structure */
   void* mumps_ptr_;
   ///@}

   /** @name Information about most recent factorization/solve */
   ///@{
   /** Number of negative eigenvalues */
   Index negevals_;
   ///@}

   /** @name Initialization flags */
   ///@{
   /** Flag indicating if internal data is initialized.
    *  For initialization, this object needs to have seen a matrix.
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

   /** Percent increase in memory */
   Index mem_percent_;

   /** Permutation and scaling method in MUMPS */
   Index mumps_permuting_scaling_;

   /** Pivot order in MUMPS. */
   Index mumps_pivot_order_;

   /** Scaling in MUMPS */
   Index mumps_scaling_;

   /** Threshold in MUMPS to state that a constraint is linearly dependent */
   Number mumps_dep_tol_;

   /** Flag indicating whether the TNLP with identical structure has
    *  already been solved before.
    */
   bool warm_start_same_structure_;
   ///@}

   /** Flag indicating if symbolic factorization has already been called */
   bool have_symbolic_factorization_;

   /** @name Internal functions */
   ///@{
   /** Call MUMPS (job=1) to perform symbolic manipulations, and reserve
    *  memory.
    */
   ESymSolverStatus SymbolicFactorization();

   /** Call MUMPS (job=2) to factorize the Matrix.
    *  It is assumed that the first nonzeros_ element of a_ contain the values
    *  of the matrix to be factorized.
    */
   ESymSolverStatus Factorization(
      bool  check_NegEVals,
      Index numberOfNegEVals
   );

   /** Call MUMPS (job=3) to do the solve. */
   ESymSolverStatus Solve(
      Index   nrhs,
      Number* rhs_vals
   );
   ///@}
};

} // namespace Ipopt
#endif
