// Copyright (C) 2009, 2010 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter, Sanjeeb Dash     IBM    2009-06-18
//                  (based on IpTSymLinearSolver.hpp rev 1019)

#include "IpParTSymLinearSolver.hpp"
#include "IpParTripletHelper.hpp"
#include "IpTripletHelper.hpp"
#include "IpBlas.hpp"

#include "IpMpi.hpp"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  ParTSymLinearSolver::ParTSymLinearSolver
  (SmartPtr<SparseSymLinearSolverInterface> solver_interface,
   SmartPtr<TSymScalingMethod> scaling_method,
   bool call_solverinterface_on_all_procs, /* = false */
   bool call_scalingmethod_on_all_procs /* = false */)
      :
      SymLinearSolver(),
      atag_(0),
      dim_(0),
      nonzeros_triplet_(0),
      nonzeros_compressed_(0),
      have_structure_(false),
      initialized_(false),

      solver_interface_(solver_interface),
      scaling_method_(scaling_method),
      scaling_factors_(NULL),
      airn_(NULL),
      ajcn_(NULL),
      recvcounts_(NULL),
      displs_(NULL),
      call_solverinterface_on_all_procs_(call_solverinterface_on_all_procs),
      call_scalingmethod_on_all_procs_(call_scalingmethod_on_all_procs)
  {
    DBG_START_METH("ParTSymLinearSolver::ParTSymLinearSolver()",dbg_verbosity);
    DBG_ASSERT(IsValid(solver_interface));

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc_);
  }

  ParTSymLinearSolver::~ParTSymLinearSolver()
  {
    DBG_START_METH("ParTSymLinearSolver::~ParTSymLinearSolver()",
                   dbg_verbosity);
    delete [] airn_;
    delete [] ajcn_;
    delete [] scaling_factors_;
    delete [] recvcounts_;
    delete [] displs_;
  }

  void ParTSymLinearSolver::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {}

  bool ParTSymLinearSolver::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    if (IsValid(scaling_method_)) {
      options.GetBoolValue("linear_scaling_on_demand",
                           linear_scaling_on_demand_, prefix);
    }
    else {
      linear_scaling_on_demand_ = false;
    }
    // This option is registered by OrigIpoptNLP
    options.GetBoolValue("warm_start_same_structure",
                         warm_start_same_structure_, prefix);

    int retval;
    if (HaveIpData()) {
      if (call_solverinterface_on_all_procs_) {
        retval = solver_interface_->Initialize(Jnlst(), IpNLP(), IpData(),
                                               IpCq(), options, prefix);
      }
      else {
        if (my_rank_==0) {
          retval = solver_interface_->Initialize(Jnlst(), IpNLP(), IpData(),
                                                 IpCq(), options, prefix);
        }
        MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
      }
    }
    else {
      if (call_solverinterface_on_all_procs_) {
        retval = solver_interface_->ReducedInitialize(Jnlst(), options, prefix);
      }
      else {
        if (my_rank_==0) {
          retval =
            solver_interface_->ReducedInitialize(Jnlst(), options, prefix);
        }
        MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
      }
    }
    if (!retval) {
      return false;
    }

    if (!warm_start_same_structure_) {
      // Reset all private data
      atag_=0;
      dim_=0;
      nonzeros_triplet_=0;
      nonzeros_compressed_=0;
      have_structure_=false;

      if (call_solverinterface_on_all_procs_) {
        matrix_format_ = solver_interface_->MatrixFormat();
      }
      else {
        int tmp;
        if (my_rank_==0) {
          matrix_format_ = solver_interface_->MatrixFormat();
          tmp = matrix_format_;
        }
        MPI_Bcast(&tmp, 1, MPI_INT, 0, MPI_COMM_WORLD);
        matrix_format_ = (SparseSymLinearSolverInterface::EMatrixFormat)tmp;
      }
      switch (matrix_format_) {
      case SparseSymLinearSolverInterface::CSR_Format_0_Offset:
        triplet_to_csr_converter_ = new TripletToCSRConverter(0);
        break;
      case SparseSymLinearSolverInterface::CSR_Format_1_Offset:
        triplet_to_csr_converter_ = new TripletToCSRConverter(1);
        break;
      case SparseSymLinearSolverInterface::CSR_Full_Format_0_Offset:
        triplet_to_csr_converter_ = new TripletToCSRConverter(0, false,
                                    TripletToCSRConverter::Full_Format);
        break;
      case SparseSymLinearSolverInterface::CSR_Full_Format_1_Offset:
        triplet_to_csr_converter_ = new TripletToCSRConverter(1, false,
                                    TripletToCSRConverter::Full_Format);
        break;
      case SparseSymLinearSolverInterface::Triplet_Format:
        triplet_to_csr_converter_ = NULL;
        break;
      default:
        Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                       "Unhandled matrix format %d in ParTSymLinearSolver::InitializeImpl\n", matrix_format_);
        return false;
      }
    }
    else {
      ASSERT_EXCEPTION(have_structure_, INVALID_WARMSTART,
                       "ParTSymLinearSolver called with warm_start_same_structure, but the internal structures are not initialized.");
    }

    // reset the initialize flag to make sure that InitializeStructure
    // is called for the linear solver
    initialized_=false;

    if (IsValid(scaling_method_) && !linear_scaling_on_demand_) {
      use_scaling_ = true;
    }
    else {
      use_scaling_ = false;
    }
    just_switched_on_scaling_ = false;

    if (IsValid(scaling_method_)) {
      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemScaling().Start();
        retval = scaling_method_->Initialize(Jnlst(), IpNLP(), IpData(),
                                             IpCq(), options, prefix);
        IpData().TimingStats().LinearSystemScaling().End();
      }
      else {
        retval = scaling_method_->ReducedInitialize(Jnlst(), options, prefix);
      }
    }
    return (bool)retval;
  }

  ESymSolverStatus
  ParTSymLinearSolver::MultiSolve(const SymMatrix& sym_A,
                                  std::vector<SmartPtr<const Vector> >& rhsV,
                                  std::vector<SmartPtr<Vector> >& solV,
                                  bool check_NegEVals,
                                  Index numberOfNegEVals)
  {
    DBG_START_METH("ParTSymLinearSolver::MultiSolve",dbg_verbosity);
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

    //DBG_ASSERT(nonzeros_triplet_== TripletHelper::GetNumberEntries(sym_A));

    // Check if the matrix has been changed
    DBG_PRINT((1,"atag_=%d sym_A->GetTag()=%d\n",atag_,sym_A.GetTag()));
    bool new_matrix = sym_A.HasChanged(atag_);
    atag_ = sym_A.GetTag();

    // If a new matrix is encountered, get the array for storing the
    // entries from the linear solver interface, fill in the new
    // values, compute the new scaling factors (if required), and
    // scale the matrix
    if (new_matrix || just_switched_on_scaling_) {
      GiveMatrixToSolver(true, sym_A);
      new_matrix = true;
    }

    // To make parallel vectors work, all processors (for now) need to
    // get all values

    // Retrieve the right hand sides and scale if required
    Index nrhs = (Index)rhsV.size();
    double* rhs_vals = new double[dim_*nrhs];
    for (Index irhs=0; irhs<nrhs; irhs++) {
      ParTripletHelper::FillAllValuesFromVector(dim_, *rhsV[irhs],
          &rhs_vals[irhs*(dim_)]);
      if (Jnlst().ProduceOutput(J_MOREMATRIX, J_LINEAR_ALGEBRA)) {
        Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                       "Right hand side %d in ParTSymLinearSolver:\n", irhs);
        for (Index i=0; i<dim_; i++) {
          Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                         "Trhs[%5d,%5d] = %23.16e\n", irhs, i,
                         rhs_vals[irhs*(dim_)+i]);
        }
      }
      if (my_rank_==0 && use_scaling_) {
        if (HaveIpData()) {
          IpData().TimingStats().LinearSystemScaling().Start();
        }
        for (Index i=0; i<dim_; i++) {
          rhs_vals[irhs*(dim_)+i] *= scaling_factors_[i];
        }
        if (HaveIpData()) {
          IpData().TimingStats().LinearSystemScaling().End();
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
      const Index* ia = NULL;
      const Index* ja = NULL;
      if (my_rank_==0) {
        if (matrix_format_==SparseSymLinearSolverInterface::Triplet_Format) {
          ia = airn_;
          ja = ajcn_;
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
      }
      if (call_solverinterface_on_all_procs_) {
        retval = solver_interface_->MultiSolve(new_matrix, ia, ja,
                                               nrhs, rhs_vals, check_NegEVals,
                                               numberOfNegEVals);
      }
      else {
        if (my_rank_==0) {
          retval = solver_interface_->MultiSolve(new_matrix, ia, ja,
                                                 nrhs, rhs_vals, check_NegEVals,
                                                 numberOfNegEVals);
        }
        MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&numberOfNegEVals, 1, MPI_INT, 0, MPI_COMM_WORLD);
      }

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
        if (my_rank_==0) {
          if (use_scaling_) {
            if (HaveIpData()) {
              IpData().TimingStats().LinearSystemScaling().Start();
            }
            for (Index i=0; i<dim_; i++) {
              rhs_vals[irhs*(dim_)+i] *= scaling_factors_[i];
            }
            if (HaveIpData()) {
              IpData().TimingStats().LinearSystemScaling().End();
            }
          }
        }
        if (Jnlst().ProduceOutput(J_MOREMATRIX, J_LINEAR_ALGEBRA)) {
          Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                         "Solution %d in ParTSymLinearSolver:\n", irhs);
          for (Index i=0; i<dim_; i++) {
            Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                           "Tsol[%5d,%5d] = %23.16e\n", irhs, i,
                           rhs_vals[irhs*(dim_)+i]);
          }
        }
        MPI_Bcast(&rhs_vals[irhs*dim_], dim_, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        ParTripletHelper::PutAllValuesInVector(dim_, &rhs_vals[irhs*(dim_)],
                                               *solV[irhs]);
      }
    }

    delete[] rhs_vals;

    return retval;
  }

  // Initialize the local copy of the positions of the nonzero
  // elements
  ESymSolverStatus
  ParTSymLinearSolver::InitializeStructure(const SymMatrix& sym_A)
  {
    DBG_START_METH("ParTSymLinearSolver::InitializeStructure",
                   dbg_verbosity);
    DBG_ASSERT(!initialized_);

    ESymSolverStatus retval;

    // have_structure_ is already true if this is a warm start for a
    // problem with identical structure
    if (!have_structure_) {

      dim_ = sym_A.Dim();

      local_nonzeros_triplet_ = ParTripletHelper::GetNumberEntries(sym_A);
      int itmp;
      MPI_Allreduce(&local_nonzeros_triplet_, &itmp, 1,
                    MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      nonzeros_triplet_ = itmp;

      Index* airn_local = new Index[local_nonzeros_triplet_];
      Index* ajcn_local = new Index[local_nonzeros_triplet_];
      ParTripletHelper::FillRowCol(local_nonzeros_triplet_, sym_A, airn_local, ajcn_local);

      if (my_rank_==0) {
        delete [] airn_;
        delete [] ajcn_;
        airn_ = new Index[nonzeros_triplet_];
        ajcn_ = new Index[nonzeros_triplet_];
        delete [] recvcounts_;
        delete [] displs_;
        recvcounts_ = new int[num_proc_];
        displs_ = new int[num_proc_];
      }

      MPI_Gather(&local_nonzeros_triplet_, 1, MPI_INT,
                 recvcounts_, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if (my_rank_==0) {
        displs_[0] = 0;
        for (int i=1; i<num_proc_; i++) {
          displs_[i] = displs_[i-1] + recvcounts_[i-1];
        }
      }

      MPI_Gatherv(airn_local, local_nonzeros_triplet_, MPI_INT,
                  airn_, recvcounts_, displs_, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Gatherv(ajcn_local, local_nonzeros_triplet_, MPI_INT,
                  ajcn_, recvcounts_, displs_, MPI_INT, 0, MPI_COMM_WORLD);

      delete [] airn_local;
      delete [] ajcn_local;

      const Index *ia = NULL;
      const Index *ja = NULL;
      Index nonzeros = -1;
      if (my_rank_==0) {
        // TODO: Catch return code for ALL processors

        // If the solver wants the compressed format, the converter has to
        // be initialized
        if (matrix_format_ == SparseSymLinearSolverInterface::Triplet_Format) {
          ia = airn_;
          ja = ajcn_;
          nonzeros = nonzeros_triplet_;
        }
        else {
          if (HaveIpData()) {
            IpData().TimingStats().LinearSystemStructureConverter().Start();
            IpData().TimingStats().LinearSystemStructureConverterInit().Start();
          }
          nonzeros_compressed_ =
            triplet_to_csr_converter_->InitializeConverter(dim_, nonzeros_triplet_,
                airn_, ajcn_);
          if (HaveIpData()) {
            IpData().TimingStats().LinearSystemStructureConverterInit().End();
          }
          ia = triplet_to_csr_converter_->IA();
          ja = triplet_to_csr_converter_->JA();
          if (HaveIpData()) {
            IpData().TimingStats().LinearSystemStructureConverter().End();
          }
          nonzeros = nonzeros_compressed_;
        }

      }

      if (call_solverinterface_on_all_procs_) {
        retval = solver_interface_->InitializeStructure(dim_, nonzeros, ia, ja);
      }
      else {
        if (my_rank_==0) {
          retval = solver_interface_->InitializeStructure(dim_, nonzeros, ia, ja);
        }
        MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
      }
      if (retval != SYMSOLVER_SUCCESS) {
        return retval;
      }

      if (my_rank_==0 || call_scalingmethod_on_all_procs_) {
        // Get space for the scaling factors
        delete [] scaling_factors_;
        if (IsValid(scaling_method_)) {
          if (HaveIpData()) {
            IpData().TimingStats().LinearSystemScaling().Start();
          }
          scaling_factors_ = new double[dim_];
          if (HaveIpData()) {
            IpData().TimingStats().LinearSystemScaling().End();
          }
        }
      }

      have_structure_ = true;
    }
    else {
      ASSERT_EXCEPTION(dim_==sym_A.Dim(), INVALID_WARMSTART,
                       "ParTSymLinearSolver called with warm_start_same_structure, but the problem is solved for the first time.");
      const Index *ia = NULL;
      const Index *ja = NULL;
      Index nonzeros = -1;
      if (my_rank_==0) {
        // This is a warm start for identical structure, so we don't need to
        // recompute the nonzeros location arrays
        if (matrix_format_ == SparseSymLinearSolverInterface::Triplet_Format) {
          ia = airn_;
          ja = ajcn_;
          nonzeros = nonzeros_triplet_;
        }
        else {
          IpData().TimingStats().LinearSystemStructureConverter().Start();
          ia = triplet_to_csr_converter_->IA();
          ja = triplet_to_csr_converter_->JA();
          IpData().TimingStats().LinearSystemStructureConverter().End();
          nonzeros = nonzeros_compressed_;
        }
      }
      if (call_solverinterface_on_all_procs_) {
        retval = solver_interface_->InitializeStructure(dim_, nonzeros, ia, ja);
      }
      else {
        if (my_rank_==0) {
          retval = solver_interface_->InitializeStructure(dim_, nonzeros, ia, ja);
        }
        MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
      }
    }
    initialized_=true;
    return retval;
  }

  Index ParTSymLinearSolver::NumberOfNegEVals() const
  {
    DBG_START_METH("ParTSymLinearSolver::NumberOfNegEVals",dbg_verbosity);

    int retval;
    if (call_solverinterface_on_all_procs_) {
      retval = solver_interface_->NumberOfNegEVals();
    }
    else {
      if (my_rank_==0) {
        retval = solver_interface_->NumberOfNegEVals();
      }
      MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    return retval;
  }

  bool ParTSymLinearSolver::IncreaseQuality()
  {
    DBG_START_METH("ParTSymLinearSolver::IncreaseQuality",dbg_verbosity);

    int retval;

    if (IsValid(scaling_method_) && !use_scaling_ &&
        linear_scaling_on_demand_) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Switching on scaling of the linear system (on demand).\n");
      IpData().Append_info_string("Mc");
      use_scaling_ = true;
      just_switched_on_scaling_ = true;
      retval = 1;
//    TODO...
//    int itmp = just_switched_on_scaling_;
//    MPI_Bcast(&itmp, 1, MPI_INT, 0, MPI_COMM_WORLD);
//    just_switched_on_scaling_ = itmp;
    }
    else {
      if (call_solverinterface_on_all_procs_) {
        retval = solver_interface_->IncreaseQuality();
      }
      else {
        if (my_rank_==0) {
          retval = solver_interface_->IncreaseQuality();
        }
        MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
      }
    }
    return retval;
  }

  bool ParTSymLinearSolver::ProvidesInertia() const
  {
    DBG_START_METH("ParTSymLinearSolver::ProvidesInertia",dbg_verbosity);

    int retval;
    if (call_solverinterface_on_all_procs_) {
      retval = solver_interface_->ProvidesInertia();
    }
    else {
      if (my_rank_==0) {
        retval = solver_interface_->ProvidesInertia();
      }
      MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    return retval;
  }

  void ParTSymLinearSolver::GiveMatrixToSolver(bool new_matrix,
      const SymMatrix& sym_A)
  {
    DBG_START_METH("ParTSymLinearSolver::GiveMatrixToSolver",dbg_verbosity);
    DBG_PRINT((1,"new_matrix = %d\n",new_matrix));

    double* pa = NULL;
    double* atriplet = NULL;

    if (my_rank_==0) {
      pa = solver_interface_->GetValuesArrayPtr();
      if (matrix_format_!=SparseSymLinearSolverInterface::Triplet_Format) {
        atriplet = new double[nonzeros_triplet_];
      }
      else {
        atriplet = pa;
      }
    }

    Number* local_atriplet = new Number[local_nonzeros_triplet_];
    ParTripletHelper::FillValues(local_nonzeros_triplet_, sym_A, local_atriplet);
    MPI_Gatherv(local_atriplet, local_nonzeros_triplet_, MPI_DOUBLE,
                atriplet, recvcounts_, displs_, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    delete [] local_atriplet;

// todo SCATTERED scaling
    if (use_scaling_) {
      if (my_rank_==0 || call_scalingmethod_on_all_procs_) {
        IpData().TimingStats().LinearSystemScaling().Start();
        DBG_ASSERT(scaling_factors_);
        if (new_matrix || just_switched_on_scaling_) {
          // only compute scaling factors if the matrix has not been
          // changed since the last call to this method
          bool retval =
            scaling_method_->ComputeSymTScalingFactors(dim_, nonzeros_triplet_,
                airn_, ajcn_, atriplet, scaling_factors_);
          if (!retval) {
            Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                           "Error during computation of scaling factors.\n");
            THROW_EXCEPTION(ERROR_IN_LINEAR_SCALING_METHOD, "scaling_method_->ComputeSymTScalingFactors returned false.")
          }
          // complain if not in debug mode
          if (my_rank_==0 &&
              Jnlst().ProduceOutput(J_MOREVECTOR, J_LINEAR_ALGEBRA)) {
            for (Index i=0; i<dim_; i++) {
              Jnlst().Printf(J_MOREVECTOR, J_LINEAR_ALGEBRA,
                             "scaling factor[%6d] = %22.17e\n",
                             i, scaling_factors_[i]);
            }
          }
          just_switched_on_scaling_ = false;
        }
        if (my_rank_==0) {
          for (Index i=0; i<nonzeros_triplet_; i++) {
            atriplet[i] *=
              scaling_factors_[airn_[i]-1] * scaling_factors_[ajcn_[i]-1];
          }
        }
        IpData().TimingStats().LinearSystemScaling().End();
      }
    }

    if (my_rank_==0) {

      if (matrix_format_!=SparseSymLinearSolverInterface::Triplet_Format) {
        IpData().TimingStats().LinearSystemStructureConverter().Start();
        triplet_to_csr_converter_->ConvertValues(nonzeros_triplet_, atriplet,
            nonzeros_compressed_, pa);
        IpData().TimingStats().LinearSystemStructureConverter().End();
        delete[] atriplet;
      }
    }

    int itmp = just_switched_on_scaling_;
    MPI_Bcast(&itmp, 1, MPI_INT, 0, MPI_COMM_WORLD);
    just_switched_on_scaling_ = itmp;

  }

} // namespace Ipopt
