// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Morten Lysgaard               2023-09-19
//               original version (based on MA57TSolverInterface.hpp)

#ifndef __IPLEOPARDSOLVERINTERFACE_HPP__
#define __IPLEOPARDSOLVERINTERFACE_HPP__

#include "IpSparseSymLinearSolverInterface.hpp"
#include "IpLibraryLoader.hpp"
#include "IpTypes.h"

#include <leopard.h>

namespace Ipopt {
/** Interface to the linear symmetric indefinite sparse solver Leopard, derived from
 *  SparseSymLinearSolverInterface.
 */
class LeopardSolverInterface : public SparseSymLinearSolverInterface {
 public:
  /** @name Constructor/Destructor */
  ///@{
  /** Constructor */
  LeopardSolverInterface(
  );

  /** Destructor */
  ~LeopardSolverInterface() override;
  ///@}

  bool InitializeImpl(
      const OptionsList &options,
      const std::string &prefix
  ) override;

  /** @name Methods for requesting solution of the linear system. */
  ///@{
  ESymSolverStatus InitializeStructure(
      Index dim,
      Index nonzeros,
      const Index *airn,
      const Index *ajcn
  ) override;

  Number *GetValuesArrayPtr() override;

  ESymSolverStatus MultiSolve(
      bool new_matrix,
      const Index *airn,
      const Index *ajcn,
      Index nrhs,
      Number *rhs_vals,
      bool check_NegEVals,
      Index numberOfNegEVals
  ) override;

  Index NumberOfNegEVals() const override;
  ///@}

  //* @name Options of Linear solver */
  ///@{
  bool IncreaseQuality() override;

  bool ProvidesInertia() const override {
    return true;
  }

  EMatrixFormat MatrixFormat() const override {
    return CSR_Format_0_Offset;
  }
  ///@}

  ///@{
  static void RegisterOptions(
      SmartPtr<RegisteredOptions> roptions
  );
  ///@}

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
  LeopardSolverInterface(
      const LeopardSolverInterface &
  ) = delete;

  /** Default Assignment Operator */
  void operator=(
      const LeopardSolverInterface &
  ) = delete;
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
  EigenvalueSignCount eigenvalue_sign_count_;
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
  /** IJ */
  IjResult ij_ = {.code=LeopardReturnCode_Ok, .object=nullptr};
  /** Ordering */
  OrderingResult ordering_ = {.code=LeopardReturnCode_Ok, .object=nullptr};
  /** OrderedIj */
  OrderedIjResult ordered_ij_ = {.code=LeopardReturnCode_Ok, .object=nullptr};
  /** AssemblyTree */
  AssemblyTreeResult assembly_tree_ = {.code=LeopardReturnCode_Ok, .object=nullptr};
  /** LDLFactorization */
#ifdef IPOPT_SINGLE
  LDLFactorizationResultF32 ldl_factorization_ = {.code=LeopardReturnCode_Ok, .object=nullptr};
#else
  LDLFactorizationResultF64 ldl_factorization_ = {.code=LeopardReturnCode_Ok, .object=nullptr};
#endif

  /** Pivot tolerance */
  Number pivtol_;
  /** Maximal pivot tolerance */
  Number pivtolmax_;
  /** Array ptr containing IPopt natural Number type*/
  Number* a_;
  ///@}

  /** @name Internal functions */
  ///@{
  /** TODO Call MA57AX and reserve memory for MA57 data.
   *
   *  TODO Reserve memory for iw_ and ikeep_, call MA57AD to perform
   *  symbolic manipulations, and reserve all the remaining data memory
   */
  ESymSolverStatus SymbolicFactorization(
      const Index *airn,
      const Index *ajcn
  );

  /** TODO Call MA57BX to factorize the Matrix.
   *
   *  It is assumed that the first nonzeros_ element of a_ contain the values
   *  of the matrix to be factorized.
   */
  ESymSolverStatus Factorization(
      const Index *airn,
      const Index *ajcn,
      bool check_NegEVals,
      Index numberOfNegEVals
  );

  /** Call TODO to do the backsolve. */
  ESymSolverStatus Backsolve(
      Index nrhs,
      Number *rhs_vals
  );
  ///@}
};

} // namespace Ipopt
#endif
