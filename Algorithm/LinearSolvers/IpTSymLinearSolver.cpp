// Copyright (C) 2004, 2005 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-03-17

#include "IpTSymLinearSolver.hpp"
#include "IpTripletHelper.hpp"

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

  TSymLinearSolver::TSymLinearSolver
  (SmartPtr<SparseSymLinearSolverInterface> solver_interface,
   SmartPtr<TSymScalingMethod> scaling_method)
      :
      SymLinearSolver(),
      atag_(0),
      dim_(0),
      nonzeros_triplet_(0),
      nonzeros_compressed_(0),
      initialized_(false),

      solver_interface_(solver_interface),
      scaling_method_(scaling_method),
      scaling_factors_(NULL),
      airn_(NULL),
      ajcn_(NULL)
  {
    DBG_START_METH("TSymLinearSolver::TSymLinearSolver()",dbg_verbosity);
    DBG_ASSERT(IsValid(solver_interface));
  }

  TSymLinearSolver::~TSymLinearSolver()
  {
    DBG_START_METH("TSymLinearSolver::~TSymLinearSolver()",
                   dbg_verbosity);
    delete [] airn_;
    delete [] ajcn_;
    delete [] scaling_factors_;
  }

  bool TSymLinearSolver::InitializeImpl(const OptionsList& options,
                                        const std::string& prefix)
  {
    // Reset all private data
    atag_=0;
    dim_=0;
    nonzeros_triplet_=0;
    nonzeros_compressed_=0;
    initialized_=false;

    if (!solver_interface_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                       options, prefix)) {
      return false;
    }

    matrix_format_ = solver_interface_->MatrixFormat();
    switch (matrix_format_) {
      case SparseSymLinearSolverInterface::CSR_Format_0_Offset:
      triplet_to_csr_converter_ = new TripletToCSRConverter(0);
      break;
      case SparseSymLinearSolverInterface::CSR_Format_1_Offset:
      triplet_to_csr_converter_ = new TripletToCSRConverter(1);
      break;
      case SparseSymLinearSolverInterface::Triplet_Format:
      triplet_to_csr_converter_ = NULL;
      break;
      default:
      DBG_ASSERT(false && "Invalid MatrixFormat returned from solver interface.");
      return false;
    }

    bool retval = true;
    if (IsValid(scaling_method_)) {
      retval = scaling_method_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                           options, prefix);
    }
    return retval;
  }

  ESymSolverStatus
  TSymLinearSolver::MultiSolve(const SymMatrix& sym_A,
                               std::vector<const Vector*>& rhsV,
                               std::vector<Vector*>& solV,
                               bool check_NegEVals,
                               Index numberOfNegEVals)
  {
    DBG_START_METH("TSymLinearSolver::MultiSolve",dbg_verbosity);
    DBG_ASSERT(!check_NegEVals || ProvidesInertia());

    // Check if this object has ever seen a matrix If not,
    // allocate memory of the matrix structure and copy the nonzeros
    // structure (it is assumed that this will never change).
    if (!initialized_) {
      ESymSolverStatus retval = InitializeStructure(sym_A);
      if (retval != SYMSOLVER_SUCCESS) {
        return retval;
      }
    }

    DBG_ASSERT(nonzeros_triplet_== TripletHelper::GetNumberEntries(sym_A));

    // Check if the matrix has been changed
    DBG_PRINT((1,"atag_=%d sym_A->GetTag()=%d\n",atag_,sym_A.GetTag()));
    bool new_matrix = sym_A.HasChanged(atag_);
    atag_ = sym_A.GetTag();

    // If a new matrix is encountered, get the array for storing the
    // entries from the linear solver interface, fill in the new
    // values, compute the new scaling factors (if required), and
    // scale the matrix
    if (new_matrix) {
      GiveMatrixToSolver(true, sym_A);
    }

    // Retrieve the right hand sides and scale if required
    Index nrhs = (Index)rhsV.size();
    double* rhs_vals = new double[dim_*nrhs];
    for (Index irhs=0; irhs<nrhs; irhs++) {
      TripletHelper::FillValuesFromVector(dim_, *rhsV[irhs],
                                          &rhs_vals[irhs*(dim_)]);
      if (IsValid(scaling_method_)) {
        for (Index i=0; i<dim_; i++) {
          rhs_vals[irhs*(dim_)+i] *= scaling_factors_[i];
        }
      }
    }

    bool done = false;
    // Call the linear solver through the interface to solve the
    // system.  This is repeated, if the return values is S_CALL_AGAIN
    // after the values have been restored (this might be necessary
    // for MA27 if the size of the work space arrays was not large
    // enough).
    ESymSolverStatus retval;
    while (!done) {
      const Index* ia;
      const Index* ja;
      if (matrix_format_==SparseSymLinearSolverInterface::Triplet_Format) {
        ia = airn_;
        ja = ajcn_;
      }
      else {
        ia = triplet_to_csr_converter_->IA();
        ja = triplet_to_csr_converter_->JA();
      }

      retval = solver_interface_->MultiSolve(new_matrix, ia, ja,
                                             nrhs, rhs_vals, check_NegEVals,
                                             numberOfNegEVals);
      if (retval==SYMSOLVER_CALL_AGAIN) {
        DBG_PRINT((1, "Solver interface asks to be called again.\n"));
        GiveMatrixToSolver(false, sym_A);
      }
      else {
        done = true;
      }
    }

    // If the solve was successful, unscale the solution (if required)
    // and transfer the result into the Vectors
    if (retval==SYMSOLVER_SUCCESS) {
      for (Index irhs=0; irhs<nrhs; irhs++) {
        if (IsValid(scaling_method_)) {
          for (Index i=0; i<dim_; i++) {
            rhs_vals[irhs*(dim_)+i] *= scaling_factors_[i];
          }
        }
        TripletHelper::PutValuesInVector(dim_, &rhs_vals[irhs*(dim_)],
                                         *solV[irhs]);
      }
    }

    delete[] rhs_vals;

    return retval;
  }

  // Initialize the local copy of the positions of the nonzero
  // elements
  ESymSolverStatus
  TSymLinearSolver::InitializeStructure(const SymMatrix& sym_A)
  {
    DBG_START_METH("TSymLinearSolver::InitializeStructure",
                   dbg_verbosity);
    DBG_ASSERT(!initialized_);

    dim_ = sym_A.Dim();
    nonzeros_triplet_ = TripletHelper::GetNumberEntries(sym_A);

    delete [] airn_;
    delete [] ajcn_;
    airn_ = new Index[nonzeros_triplet_];
    ajcn_ = new Index[nonzeros_triplet_];

    TripletHelper::FillRowCol(nonzeros_triplet_, sym_A, airn_, ajcn_);

    // If the solver wants the compressed format, the converter has to
    // be initialized
    const Index *ia;
    const Index *ja;
    Index nonzeros;
    if (matrix_format_ == SparseSymLinearSolverInterface::Triplet_Format) {
      ia = airn_;
      ja = ajcn_;
      nonzeros = nonzeros_triplet_;
    }
    else {
      nonzeros_compressed_ =
        triplet_to_csr_converter_->InitializeConverter(dim_, nonzeros_triplet_,
            airn_, ajcn_);
      ia = triplet_to_csr_converter_->IA();
      ja = triplet_to_csr_converter_->JA();
      nonzeros = nonzeros_compressed_;
    }

    ESymSolverStatus retval =
      solver_interface_->InitializeStructure(dim_, nonzeros, ia, ja);
    if (retval != SYMSOLVER_SUCCESS) {
      return retval;
    }

    // Get space for the scaling factors
    delete [] scaling_factors_;
    if (IsValid(scaling_method_)) {
      scaling_factors_ = new double[dim_];
    }

    initialized_ = true;
    return retval;
  }

  Index TSymLinearSolver::NumberOfNegEVals() const
  {
    DBG_START_METH("TSymLinearSolver::NumberOfNegEVals",dbg_verbosity);
    return solver_interface_->NumberOfNegEVals();
  }

  bool TSymLinearSolver::IncreaseQuality()
  {
    DBG_START_METH("TSymLinearSolver::IncreaseQuality",dbg_verbosity);

    return solver_interface_->IncreaseQuality();
  }

  bool TSymLinearSolver::ProvidesInertia() const
  {
    DBG_START_METH("TSymLinearSolver::ProvidesInertia",dbg_verbosity);

    return solver_interface_->ProvidesInertia();
  }

  void TSymLinearSolver::GiveMatrixToSolver(bool new_matrix,
      const SymMatrix& sym_A)
  {
    DBG_START_METH("TSymLinearSolver::GiveMatrixToSolver",dbg_verbosity);
    DBG_PRINT((1,"new_matrix = %d\n",new_matrix));

    double* pa = solver_interface_->GetValuesArrayPtr();
    double* atriplet;

    if (matrix_format_!=SparseSymLinearSolverInterface::Triplet_Format) {
      atriplet = new double[nonzeros_triplet_];
    }
    else {
      atriplet = pa;
    }

    //DBG_PRINT_MATRIX(3, "Aunscaled", sym_A);
    TripletHelper::FillValues(nonzeros_triplet_, sym_A, atriplet);
    if (DBG_VERBOSITY()>=3) {
      for (Index i=0; i<nonzeros_triplet_; i++) {
        DBG_PRINT((3, "KKTunscaled(%6d,%6d) = %24.16e\n", airn_[i], ajcn_[i], atriplet[i]));
      }
    }

    if (IsValid(scaling_method_)) {
      DBG_ASSERT(scaling_factors_);
      if (new_matrix) {
        // only compute scaling factors if the matrix has not been
        // changed since the last call to this method
        bool retval =
          scaling_method_->ComputeSymTScalingFactors(dim_, nonzeros_triplet_,
              airn_, ajcn_,
              atriplet, scaling_factors_);
        DBG_ASSERT(retval);
        retval = false; // is added to make sure compiles doesn't
        // complain if not in debug mode
        if (Jnlst().ProduceOutput(J_MOREVECTOR, J_LINEAR_ALGEBRA)) {
          for (Index i=0; i<dim_; i++) {
            Jnlst().Printf(J_MOREVECTOR, J_LINEAR_ALGEBRA,
                           "scaling factor[%6d] = %22.17e\n",
                           i, scaling_factors_[i]);
          }
        }
      }
      for (Index i=0; i<nonzeros_triplet_; i++) {
        atriplet[i] *=
          scaling_factors_[airn_[i]-1] * scaling_factors_[ajcn_[i]-1];
      }
      if (DBG_VERBOSITY()>=3) {
        for (Index i=0; i<nonzeros_triplet_; i++) {
          DBG_PRINT((3, "KKTscaled(%6d,%6d) = %24.16e\n", airn_[i], ajcn_[i], atriplet[i]));
        }
      }
    }

    if (matrix_format_!=SparseSymLinearSolverInterface::Triplet_Format) {
      triplet_to_csr_converter_->ConvertValues(nonzeros_triplet_, atriplet,
          nonzeros_compressed_, pa);
      delete[] atriplet;
    }

  }

} // namespace Ipopt
