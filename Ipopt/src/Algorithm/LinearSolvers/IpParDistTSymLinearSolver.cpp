// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter, Sanjeeb Dash     IBM    2009-07-14
//                  (based on IpTSymLinearSolver.hpp rev 1019)

#include "IpParDistTSymLinearSolver.hpp"
#include "IpParTripletHelper.hpp"
#include "IpTripletHelper.hpp"
#include "IpBlas.hpp"

#include "IpMpi.hpp"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  ParDistTSymLinearSolver::ParDistTSymLinearSolver
  (SmartPtr<SparseSymLinearSolverInterface> solver_interface)
      :
      SymLinearSolver(),
      atag_(0),
      local_dim_(-1),
      have_structure_(false),
      initialized_(false),

      solver_interface_(solver_interface),
      airn_local_(NULL),
      ajcn_local_(NULL)
  {
    DBG_START_METH("ParDistTSymLinearSolver::ParDistTSymLinearSolver()",dbg_verbosity);
    DBG_ASSERT(IsValid(solver_interface));

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc_);
  }

  ParDistTSymLinearSolver::~ParDistTSymLinearSolver()
  {
    DBG_START_METH("ParDistTSymLinearSolver::~ParDistTSymLinearSolver()",
                   dbg_verbosity);
    delete [] airn_local_;
    delete [] ajcn_local_;
  }

  void ParDistTSymLinearSolver::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {}

  bool ParDistTSymLinearSolver::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    // This option is registered by OrigIpoptNLP
    options.GetBoolValue("warm_start_same_structure",
                         warm_start_same_structure_, prefix);

    bool retval;
    if (HaveIpData()) {
      retval = solver_interface_->Initialize(Jnlst(), IpNLP(), IpData(),
                                             IpCq(), options, prefix);
    }
    else {
      retval = solver_interface_->ReducedInitialize(Jnlst(), options, prefix);
    }
    if (!retval) {
      return false;
    }

    if (!warm_start_same_structure_) {
      // Reset all private data
      atag_=0;
      dim_=0;
      nonzeros_triplet_local_=0;
      nonzeros_compressed_local_=0;
      have_structure_=false;

      matrix_format_ = solver_interface_->MatrixFormat();
      switch (matrix_format_) {
      case SparseSymLinearSolverInterface::CSC_Format_0_Offset:
        triplet_to_csr_converter_ = new TripletToCSRConverter(0, true);
        break;
      case SparseSymLinearSolverInterface::CSC_Format_1_Offset:
        triplet_to_csr_converter_ = new TripletToCSRConverter(1, true);
        break;
      case SparseSymLinearSolverInterface::Triplet_Format:
        triplet_to_csr_converter_ = NULL;
        break;
      default:
        DBG_ASSERT(false && "Invalid MatrixFormat returned from solver interface.");
        return false;
      }
    }
    else {
      ASSERT_EXCEPTION(have_structure_, INVALID_WARMSTART,
                       "ParDistTSymLinearSolver called with warm_start_same_structure, but the internal structures are not initialized.");
    }

    // reset the initialize flag to make sure that InitializeStructure
    // is called for the linear solver
    initialized_=false;

    return retval;
  }

  ESymSolverStatus
  ParDistTSymLinearSolver::MultiSolve(const SymMatrix& sym_A,
                                      std::vector<SmartPtr<const Vector> >& rhsV,
                                      std::vector<SmartPtr<Vector> >& solV,
                                      bool check_NegEVals,
                                      Index numberOfNegEVals)
  {
    DBG_START_METH("ParDistTSymLinearSolver::MultiSolve",dbg_verbosity);
    DBG_ASSERT(!check_NegEVals || ProvidesInertia());

    // Check if this object has ever seen a matrix If not,
    // allocate memory of the matrix structure and copy the nonzeros
    // structure (it is assumed that this will never change).
    if (!initialized_) {
      ESymSolverStatus retval = InitializeStructure(sym_A, *rhsV[0]);
      if (retval != SYMSOLVER_SUCCESS) {
        return retval;
      }
    }

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
      new_matrix = true;
    }

    // To make parallel vectors work, all processors (for now) need to
    // get all values

    // Retrieve the right hand sides and scale if required
    Index nrhs = (Index)rhsV.size();
    double* rhs_vals;
    if (local_dim_==-1) {
      // Solver wants entire right hand side
      rhs_vals = new double[dim_*nrhs];
      for (Index irhs=0; irhs<nrhs; irhs++) {
        ParTripletHelper::FillAllValuesFromVector(dim_, *rhsV[irhs],
            &rhs_vals[irhs*(dim_)]);
        if (Jnlst().ProduceOutput(J_MOREMATRIX, J_LINEAR_ALGEBRA)) {
          Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                         "Right hand side %d in ParDistTSymLinearSolver:\n", irhs);
          for (Index i=0; i<dim_; i++) {
            Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                           "Trhs[%5d,%5d] = %23.16e\n", irhs, i,
                           rhs_vals[irhs*(dim_)+i]);
          }
        }
      }
    }
    else {
      rhs_vals = new double[local_dim_*nrhs];
      for (Index irhs=0; irhs<nrhs; irhs++) {
        ParTripletHelper::FillLocalValuesFromVector(local_dim_, *rhsV[irhs],
            &rhs_vals[irhs*(local_dim_)]);
        if (Jnlst().ProduceOutput(J_MOREMATRIX, J_LINEAR_ALGEBRA)) {
          Jnlst().StartDistributedOutput();
          Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                         "Right hand side %d in ParDistTSymLinearSolver on process %d:\n", irhs, my_rank_);
          for (Index i=0; i<local_dim_; i++) {
            Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                           "Trhs[%5d,%5d] = %23.16e\n", irhs, i,
                           rhs_vals[irhs*(local_dim_)+i]);
          }
          Jnlst().FinishDistributedOutput();
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
        ia = airn_local_;
        ja = ajcn_local_;
      }
      else {
        if (HaveIpData()) {
          IpData().TimingStats().LinearSystemStructureConverter().Start();
        }
        ia = triplet_to_csr_converter_->IA();
        ja = triplet_to_csr_converter_->JA();
        if (HaveIpData()) {
          IpData().TimingStats().LinearSystemStructureConverter().End();
        }
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
      if (local_dim_ == -1) {
        for (Index irhs=0; irhs<nrhs; irhs++) {
          if (Jnlst().ProduceOutput(J_MOREMATRIX, J_LINEAR_ALGEBRA)) {
            Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                           "Solution %d in ParDistTSymLinearSolver:\n", irhs);
            for (Index i=0; i<dim_; i++) {
              Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                             "Tsol[%5d,%5d] = %23.16e\n", irhs, i,
                             rhs_vals[irhs*(dim_)+i]);
            }
          }
          ParTripletHelper::PutAllValuesInVector(dim_, &rhs_vals[irhs*(dim_)],
                                                 *solV[irhs]);
        }
      }
      else {
        for (Index irhs=0; irhs<nrhs; irhs++) {
          if (Jnlst().ProduceOutput(J_MOREMATRIX, J_LINEAR_ALGEBRA)) {
            Jnlst().StartDistributedOutput();
            Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                           "Solution %d in ParDistTSymLinearSolver in process %d:\n", irhs, my_rank_);
            for (Index i=0; i<local_dim_; i++) {
              Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                             "Tsol[%5d,%5d] = %23.16e\n", irhs, i,
                             rhs_vals[irhs*(local_dim_)+i]);
            }
            Jnlst().FinishDistributedOutput();
          }
          ParTripletHelper::PutLocalValuesInVector(local_dim_, &rhs_vals[irhs*(local_dim_)],
              *solV[irhs]);
        }
      }
    }

    delete[] rhs_vals;

    return retval;
  }

  // Initialize the local copy of the positions of the nonzero
  // elements
  ESymSolverStatus
  ParDistTSymLinearSolver::InitializeStructure(const SymMatrix& sym_A,
      const Vector& sample_rhs)
  {
    DBG_START_METH("ParDistTSymLinearSolver::InitializeStructure",
                   dbg_verbosity);
    DBG_ASSERT(!initialized_);

    ESymSolverStatus retval;

    // have_structure_ is already true if this is a warm start for a
    // problem with identical structure
    if (!have_structure_) {

      dim_ = sym_A.Dim();

      nonzeros_triplet_local_ = ParTripletHelper::GetNumberEntries(sym_A);

      airn_local_ = new Index[nonzeros_triplet_local_];
      ajcn_local_ = new Index[nonzeros_triplet_local_];
      ParTripletHelper::FillRowCol(nonzeros_triplet_local_, sym_A, airn_local_, ajcn_local_);

      // If the solver wants the compressed format, the converter has to
      // be initialized
      const Index *ia_local;
      const Index *ja_local;
      Index nonzeros_local;
      if (matrix_format_ == SparseSymLinearSolverInterface::Triplet_Format) {
        ia_local = airn_local_;
        ja_local = ajcn_local_;
        nonzeros_local = nonzeros_triplet_local_;
      }
      else {
        if (HaveIpData()) {
          IpData().TimingStats().LinearSystemStructureConverter().Start();
          IpData().TimingStats().LinearSystemStructureConverterInit().Start();
        }
        // tell solver about global indices
        local_dim_ = ParTripletHelper::GetLocalNumberEntries(sample_rhs);
        Index* global_pos = new Index[local_dim_];
        Index offset = 0;
        if (matrix_format_ == SparseSymLinearSolverInterface::CSC_Format_1_Offset) {
          offset = 1;
        }
        ParTripletHelper::GetGlobalPos(local_dim_, sample_rhs, global_pos, offset);
        if (Jnlst().ProduceOutput(J_MOREMATRIX, J_LINEAR_ALGEBRA)) {
          Jnlst().StartDistributedOutput();
          Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                         "On Process %d, we have %d local rows with the following global positions:\n", my_rank_, local_dim_);
          for (Index i=0; i<local_dim_; ++i) {
            Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                           "  global_pos[%5d] = %5d\n", i, global_pos[i]);
          }
          Jnlst().FinishDistributedOutput();
        }
        bool retv = solver_interface_->SetGlobalPos(local_dim_, global_pos);
        delete [] global_pos;
        if (!retv) {
          Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                         "Selected linear solver rejected call to SetGlobalPos:\n");
          return SYMSOLVER_FATAL_ERROR;
        }

        if (Jnlst().ProduceOutput(J_MOREMATRIX, J_LINEAR_ALGEBRA)) {
          Jnlst().StartDistributedOutput();
          Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA, "On Process %d, TripletConverter gets matrix:\n");
          for (Index i=0; i<nonzeros_triplet_local_; i++) {
            Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                           "  airn[%5d] = %5d acjn[%5d] = %5d\n",
                           i, airn_local_[i], i, ajcn_local_[i]);
          }
          Jnlst().FinishDistributedOutput();
        }

        nonzeros_compressed_local_ =
          triplet_to_csr_converter_->InitializeConverter(dim_, nonzeros_triplet_local_,
              airn_local_, ajcn_local_);
        if (Jnlst().ProduceOutput(J_MOREMATRIX, J_LINEAR_ALGEBRA)) {
          ia_local = triplet_to_csr_converter_->IA();
          ja_local = triplet_to_csr_converter_->JA();
          if (matrix_format_ ==  SparseSymLinearSolverInterface::CSC_Format_1_Offset) {
            ja_local--;
          }
          Jnlst().StartDistributedOutput();
          Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA, "On Process %d, TripletConverter returns (uncompressed) structure:\n", my_rank_);
          for (Index i=0; i<dim_; i++) {
            Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                           "  ia_local[%5d] = %5d\n", i, ia_local[i]);
            for (Index j=ia_local[i]; j<ia_local[i+1]; j++) {
              Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                             "    ja_local[%5d] = %5d\n", j, ja_local[j]);
            }
          }
          Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                         "  ia_local[%5d] = %5d\n", dim_, ia_local[dim_]);
          Jnlst().FinishDistributedOutput();
        }
        Index retval = triplet_to_csr_converter_->DeleteZeroRows();
        DBG_ASSERT(retval == local_dim_);
        if (Jnlst().ProduceOutput(J_MOREMATRIX, J_LINEAR_ALGEBRA)) {
          ia_local = triplet_to_csr_converter_->IA();
          ja_local = triplet_to_csr_converter_->JA();
          if (matrix_format_ ==  SparseSymLinearSolverInterface::CSC_Format_1_Offset) {
            ja_local--;
          }
          Jnlst().StartDistributedOutput();
          Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA, "On Process %d, TripletConverter returns (compressed) structure:\n", my_rank_);
          for (Index i=0; i<local_dim_; i++) {
            Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                           "  ia_local[%5d] = %5d\n", i, ia_local[i]);
            for (Index j=ia_local[i]; j<ia_local[i+1]; j++) {
              Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                             "    ja_local[%5d] = %5d\n", j, ja_local[j]);
            }
          }
          Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                         "  ia_local[%5d] = %5d\n", local_dim_, ia_local[local_dim_]);
          Jnlst().FinishDistributedOutput();
        }
        ia_local = triplet_to_csr_converter_->IA();
        ja_local = triplet_to_csr_converter_->JA();

        if (HaveIpData()) {
          IpData().TimingStats().LinearSystemStructureConverterInit().End();
        }
        if (HaveIpData()) {
          IpData().TimingStats().LinearSystemStructureConverter().End();
        }
        nonzeros_local = nonzeros_compressed_local_;
      }

      retval = solver_interface_->InitializeStructure(dim_, nonzeros_local, ia_local, ja_local);
      if (retval != SYMSOLVER_SUCCESS) {
        return retval;
      }

      have_structure_ = true;
    }
    else {
      ASSERT_EXCEPTION(dim_==sym_A.Dim(), INVALID_WARMSTART,
                       "ParDistTSymLinearSolver called with warm_start_same_structure, but the problem is solved for the first time.");
      // This is a warm start for identical structure, so we don't need to
      // recompute the nonzeros location arrays
      const Index *ia_local;
      const Index *ja_local;
      Index nonzeros_local;
      if (matrix_format_ == SparseSymLinearSolverInterface::Triplet_Format) {
        ia_local = airn_local_;
        ja_local = ajcn_local_;
        nonzeros_local = nonzeros_triplet_local_;
      }
      else {
        IpData().TimingStats().LinearSystemStructureConverter().Start();
        ia_local = triplet_to_csr_converter_->IA();
        ja_local = triplet_to_csr_converter_->JA();
        IpData().TimingStats().LinearSystemStructureConverter().End();
        nonzeros_local = nonzeros_compressed_local_;
      }
      retval = solver_interface_->InitializeStructure(dim_, nonzeros_local, ia_local, ja_local);
    }
    initialized_=true;
    return retval;
  }

  Index ParDistTSymLinearSolver::NumberOfNegEVals() const
  {
    DBG_START_METH("ParDistTSymLinearSolver::NumberOfNegEVals",dbg_verbosity);

    return solver_interface_->NumberOfNegEVals();
  }

  bool ParDistTSymLinearSolver::IncreaseQuality()
  {
    DBG_START_METH("ParDistTSymLinearSolver::IncreaseQuality",dbg_verbosity);

    return solver_interface_->IncreaseQuality();
  }

  bool ParDistTSymLinearSolver::ProvidesInertia() const
  {
    DBG_START_METH("ParDistTSymLinearSolver::ProvidesInertia",dbg_verbosity);

    return solver_interface_->ProvidesInertia();
  }

  void ParDistTSymLinearSolver::GiveMatrixToSolver(bool new_matrix,
      const SymMatrix& sym_A)
  {
    DBG_START_METH("ParDistTSymLinearSolver::GiveMatrixToSolver",dbg_verbosity);
    DBG_PRINT((1,"new_matrix = %d\n",new_matrix));

    double* pa;
    double* atriplet_local;

    pa = solver_interface_->GetValuesArrayPtr();
    if (matrix_format_!=SparseSymLinearSolverInterface::Triplet_Format) {
      atriplet_local = new double[nonzeros_triplet_local_];
    }
    else {
      atriplet_local = pa;
    }

    ParTripletHelper::FillValues(nonzeros_triplet_local_, sym_A, atriplet_local);

    // print matrix values...
    if (Jnlst().ProduceOutput(J_MOREMATRIX, J_LINEAR_ALGEBRA)) {
      Jnlst().StartDistributedOutput();
      Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                     "Matrix values on process %d (in trilpet order):\n", my_rank_);
      for (Index i=0; i<nonzeros_triplet_local_; ++i) {
        Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                       "atrp[%7d] = %23.16e\n", i, atriplet_local[i]);
      }
      Jnlst().FinishDistributedOutput();
    }

    if (matrix_format_!=SparseSymLinearSolverInterface::Triplet_Format) {
      IpData().TimingStats().LinearSystemStructureConverter().Start();
      triplet_to_csr_converter_->ConvertValues(nonzeros_triplet_local_, atriplet_local,
          nonzeros_compressed_local_, pa);
      IpData().TimingStats().LinearSystemStructureConverter().End();
      delete[] atriplet_local;
    }

  }

} // namespace Ipopt
