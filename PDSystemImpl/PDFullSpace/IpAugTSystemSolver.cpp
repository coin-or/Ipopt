// Copyright (C) 2004, International BusinDess Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpAugTSystemSolver.hpp"
#include "IpTripletHelper.hpp"
#include "IpBlas.hpp"
#include "IpDebug.hpp"

namespace Ipopt
{

  static const Index dbg_verbosity = 0;

  AugTSystemSolver::AugTSystemSolver(SymLinearSolver& linSolver)
      :
      AugSystemSolver(),
      linsolver_(&linSolver),
      w_tag_(0),
      d_x_tag_(0),
      delta_x_(0.),
      d_s_tag_(0),
      delta_s_(0.),
      j_c_tag_(0),
      d_c_tag_(0),
      delta_c_(0.),
      j_d_tag_(0),
      d_d_tag_(0),
      delta_d_(0.),
      augsys_tag_(0),
      nx_(0),
      nd_(0),
      nc_(0),
      nnz_w_(0),
      nnz_j_c_(0),
      nnz_j_d_(0),
      augsystem_(NULL),
      initialized_(false),
      have_values_(false)
  {
    DBG_START_METH("AugTSystemSolver::AugTSystemSolver()",dbg_verbosity);
  }

  AugTSystemSolver::~AugTSystemSolver()
  {
    DBG_START_METH("AugTSystemSolver::~AugTSystemSolver()",dbg_verbosity);
  }


  bool AugTSystemSolver::InitializeImpl(const OptionsList& options,
                                        const std::string& prefix)
  {
    // TODO initialize tags and flags!

    return linsolver_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                  options, prefix);
  }


  SymLinearSolver::ESolveStatus
  AugTSystemSolver::Solve(const SymMatrix* W,
                          const Vector* D_x,
                          double delta_x,
                          const Vector* D_s,
                          double delta_s,
                          const Matrix* J_c,
                          const Vector* D_c,
                          double delta_c,
                          const Matrix* J_d,
                          const Vector* D_d,
                          double delta_d,
                          const Vector& rhs_x,
                          const Vector& rhs_s,
                          const Vector& rhs_c,
                          const Vector& rhs_d,
                          Vector& sol_x,
                          Vector& sol_s,
                          Vector& sol_c,
                          Vector& sol_d,
                          bool check_NegEVals,
                          Index numberOfNegEVals)
  {
    DBG_START_METH("AugTSystemSolver::Solve",dbg_verbosity);

    // If the internal copy of the augmented system matrix has not yet been
    // initialized, do it now
    if (!initialized_) {
      InternalInitialize(W, J_c, J_d);
    }

    // Use the Set method to update (if necessary) the internal copy of
    // the augmented system matrix
    Set(W, D_x, delta_x, D_s, delta_s, J_c, D_c, delta_c, J_d, D_d, delta_d);

    // Sanity checks
    DBG_ASSERT(rhs_x.Dim()==nx_);
    DBG_ASSERT(rhs_s.Dim()==nd_);
    DBG_ASSERT(rhs_c.Dim()==nc_);
    DBG_ASSERT(rhs_d.Dim()==nd_);
    DBG_ASSERT(sol_x.Dim()==nx_);
    DBG_ASSERT(sol_s.Dim()==nd_);
    DBG_ASSERT(sol_c.Dim()==nc_);
    DBG_ASSERT(sol_d.Dim()==nd_);

    // Now construct the overall right hand side vector that will be passed
    // to the linear solver
    aug_vec_space_ = new DenseVectorSpace(nx_+nd_+nc_+nd_);
    SmartPtr<DenseVector> aug_rhs = aug_vec_space_->MakeNewDenseVector();
    aug_rhs->CopyToPos(0, rhs_x);
    aug_rhs->CopyToPos(nx_, rhs_s);
    aug_rhs->CopyToPos(nx_+nd_, rhs_c);
    aug_rhs->CopyToPos(nx_+nd_+nc_, rhs_d);

    Jnlst().PrintMatrix(J_MATRIX, J_LINEAR_ALGEBRA, "KKT", *augsystem_);
    Jnlst().PrintVector(J_MOREVECTOR, J_LINEAR_ALGEBRA, "RHS", *aug_rhs);
    // Call the linear solver
    SmartPtr<DenseVector> aug_sol = aug_vec_space_->MakeNewDenseVector();
    SymLinearSolver::ESolveStatus retval;
    retval = linsolver_->Solve(*augsystem_, *aug_rhs, *aug_sol,
                               check_NegEVals, numberOfNegEVals);
    Jnlst().PrintVector(J_MOREVECTOR, J_LINEAR_ALGEBRA, "SOL", *aug_sol);

    if (retval==SymLinearSolver::S_SUCCESS) {
      // Copy the solution into the individual sol vectors
      aug_sol->CopyFromPos(0, sol_x);
      aug_sol->CopyFromPos(nx_, sol_s);
      aug_sol->CopyFromPos(nx_+nd_, sol_c);
      aug_sol->CopyFromPos(nx_+nd_+nc_, sol_d);
    }

    return retval;
  }

  void AugTSystemSolver::InternalInitialize(const SymMatrix* W,
      const Matrix* J_c,
      const Matrix* J_d)
  {
    DBG_START_METH("AugTSystemSolver::InternalInitialize",dbg_verbosity);
    // Sanity checks
    DBG_ASSERT(!initialized_);

    // Obtain data about the size of the problem
    const GenTMatrix* tmatrix_j_c =
      dynamic_cast<const GenTMatrix*>(J_c);
    DBG_ASSERT(tmatrix_j_c);
    const GenTMatrix* tmatrix_j_d =
      dynamic_cast<const GenTMatrix*>(J_d);
    DBG_ASSERT(tmatrix_j_d);

    // Get the dimensions of the problem
    nx_ = J_c->NCols(); // Number of optimization variables
    nc_ = J_c->NRows(); // Number of equality constraints
    nd_ = J_d->NRows(); // Number of inequality constraints
    DBG_ASSERT(nx_==J_d->NCols());

    nnz_j_c_ = tmatrix_j_c->Nonzeros();
    nnz_j_d_ = tmatrix_j_d->Nonzeros();

    // If the W matrix is set to NULL, the we assume that the upper left
    // part only consists of a diagonal
    const SymTMatrix* tmatrix_w;
    if (W==NULL) {
      nnz_w_ = 0;
    }
    else {
      DBG_ASSERT(nx_==W->Dim());
      tmatrix_w = dynamic_cast<const SymTMatrix*>(W);
      DBG_ASSERT(tmatrix_w);
      nnz_w_ = tmatrix_w->Nonzeros();
    }

    // allocate the SymTMatrix for storing the augmented system
    // For this, construct arrays with the positions of the nonzeros
    Index nnz_augsys = nx_ + nd_ + nnz_w_ + nnz_j_c_ + nnz_j_d_ + nd_
                       + nc_ + nd_;
    Index* irn = new Index[nnz_augsys];
    Index* jcn = new Index[nnz_augsys];

    Index ioff = 0;
    for(Index i=0; i<nx_+nd_; i++) {
      irn[ioff+i] = i+1;
      jcn[ioff+i] = i+1;
    }
    ioff += nx_+nd_;

    if (nnz_w_>0) {
      // This is skipped if W is NULL
      const Index* irn_w = tmatrix_w->Irows();
      const Index* jcn_w = tmatrix_w->Jcols();
      for(Index i=0; i<nnz_w_; i++) {
        irn[ioff+i] = irn_w[i];
        jcn[ioff+i] = jcn_w[i];
      }
      ioff+=nnz_w_;
    }

    const Index* irows_c = tmatrix_j_c->Irows();
    const Index* jcols_c = tmatrix_j_c->Jcols();
    for(Index i=0; i<nnz_j_c_; i++) {
      irn[ioff+i] = irows_c[i]+nx_+nd_;
      jcn[ioff+i] = jcols_c[i];
    }
    ioff+=nnz_j_c_;

    const Index* irows_d = tmatrix_j_d->Irows();
    const Index* jcols_d = tmatrix_j_d->Jcols();
    for(Index i=0; i<nnz_j_d_; i++) {
      irn[ioff+i] = irows_d[i]+nx_+nd_+nc_;
      jcn[ioff+i] = jcols_d[i];
    }
    ioff+=nnz_j_d_;

    for(Index i=0; i<nd_; i++) {
      irn[ioff+i] = 1+i+nx_+nd_+nc_;
      jcn[ioff+i] = 1+i+nx_;
    }
    ioff+=nd_;

    for(Index i=0; i<nc_+nd_; i++) {
      irn[ioff+i] = nx_+nd_+i+1;
      jcn[ioff+i] = nx_+nd_+i+1;
    }
    ioff+=nc_+nd_;
    DBG_ASSERT(ioff==nnz_augsys);

    DBG_ASSERT(IsNull(augsystem_));
    aug_mat_space_ = new SymTMatrixSpace(nx_+nd_+nc_+nd_, nnz_augsys, irn, jcn);
    augsystem_ = aug_mat_space_->MakeNewSymTMatrix();
    delete [] irn;
    delete [] jcn;

    augsys_tag_ = augsystem_->GetTag();
    initialized_ = true;
  }

  void AugTSystemSolver::Set(const SymMatrix* W,
                             const Vector* D_x,
                             double delta_x,
                             const Vector* D_s,
                             double delta_s,
                             const Matrix* J_c,
                             const Vector* D_c,
                             double delta_c,
                             const Matrix* J_d,
                             const Vector* D_d,
                             double delta_d)
  {
    DBG_START_METH("AugTSystemSolver::Set",dbg_verbosity);
    // Sanity checks
    DBG_ASSERT(initialized_);
    // if W was NULL for the initialization, but is not NULL now
    DBG_ASSERT(W==NULL || nnz_w_>0);
    // Make sure noone from the outside messed with our matrix
    DBG_ASSERT(!augsystem_->HasChanged(augsys_tag_));

    // Figure out which parts of the matrix have to be (re)computed.
    bool change_in_D_x = true;
    bool change_in_D_s = true;
    bool change_in_W = true;
    bool change_in_J_c = true;
    bool change_in_J_d = true;
    bool change_in_D_c = true;
    bool change_in_D_d = true;
    if (have_values_) {
      if (delta_x_==delta_x) {
        change_in_D_x = ChangeInPart(D_x,d_x_tag_);
      }
      if (delta_s_==delta_s) {
        change_in_D_s = ChangeInPart(D_s,d_s_tag_);
      }
      change_in_W = ChangeInPart(W,w_tag_);
      change_in_J_c = ChangeInPart(J_c,j_c_tag_);
      change_in_J_d = ChangeInPart(J_d,j_d_tag_);
      if (delta_c_==delta_c) {
        change_in_D_c = ChangeInPart(D_c,d_c_tag_);
      }
      if (delta_d_==delta_d) {
        change_in_D_d = ChangeInPart(D_d,d_d_tag_);
      }
    }

    // If there has been any change, request the internal Values array from
    // the augmented system matrix (which will mark that matrix as changed)
    // and recompute the parts of that array that have changed.
    if (change_in_D_x || change_in_D_s || change_in_W || change_in_J_c ||
        change_in_J_d || change_in_D_c || change_in_D_d) {
      // The following will mark the sp_augsystem_ matrix as changed, so that
      // the linear solver object knows that it has to do the factorization
      // again
      Number* val_augsys = augsystem_->Values();

      // Let's start with the upper left corner diagonals
      Index ioff=0;
      if (change_in_D_x) {
        DBG_PRINT((2,"Computing D_x entries.\n"));
        FillDPart(nx_, delta_x, D_x, val_augsys+ioff);
      }
      ioff+=nx_;
      if (change_in_D_s) {
        DBG_PRINT((2,"Computing D_s entries.\n"));
        FillDPart(nd_, delta_s, D_s, val_augsys+ioff);
      }
      ioff+=nd_;
      // Now let's take care of the W matrix (if it is part of the structure)
      if (change_in_W) {
        DBG_PRINT((2,"Computing W entries.\n"));
        FillWPart(nnz_w_, W, val_augsys+ioff);
      }
      ioff+=nnz_w_;
      // Now the Jacobians
      if (change_in_J_c) {
        DBG_PRINT((2,"Computing J_c entries.\n"));
        FillJPart(nnz_j_c_, J_c, val_augsys+ioff);
      }
      ioff+=nnz_j_c_;
      if (change_in_J_d) {
        DBG_PRINT((2,"Computing J_d entries.\n"));
        FillJPart(nnz_j_d_, J_d, val_augsys+ioff);
      }
      ioff+=nnz_j_d_;
      // Here we write the values for the Identity (Jacobian w.r.t.
      // inequality slack variables)
      if (!have_values_) {
        DBG_PRINT((2,"Initializing -Identity entries.\n"));
        Number mone = -1.;
        IpBlasDcopy(nd_, &mone, 0, val_augsys+ioff, 1);
      }
      ioff+=nd_;
      // Finally the lower right corner diagonals
      if (change_in_D_c) {
        DBG_PRINT((2,"Computing D_c entries.\n"));
        FillDPart(nc_, -delta_c, D_c, val_augsys+ioff);
      }
      ioff+=nc_;
      if (change_in_D_d) {
        DBG_PRINT((2,"Computing D_d entries.\n"));
        FillDPart(nd_, -delta_d, D_d, val_augsys+ioff);
      }
      ioff+=nd_;
      DBG_ASSERT(ioff==augsystem_->Nonzeros());
    } //if (change_in_D_x ...

    // Store the current values and tag for the next call of Set
    delta_x_ = delta_x;
    delta_s_ = delta_s;
    delta_c_ = delta_c;
    delta_d_ = delta_d;

    d_x_tag_ = NewTag(D_x);
    d_s_tag_ = NewTag(D_s);
    w_tag_ = NewTag(W);
    j_c_tag_ = NewTag(J_c);
    j_d_tag_ = NewTag(J_d);
    d_c_tag_ = NewTag(D_c);
    d_d_tag_ = NewTag(D_d);

    have_values_ = true;
    augsys_tag_ = augsystem_->GetTag();
  }

  Index AugTSystemSolver::NumberOfNegEVals() const
  {
    DBG_ASSERT(initialized_);
    return linsolver_->NumberOfNegEVals();
  }

  bool AugTSystemSolver::ProvidesInertia() const
  {
    return linsolver_->ProvidesInertia();
  }

  bool AugTSystemSolver::ChangeInPart(const TaggedObject* Obj,
                                      const TaggedObject::Tag tag)
  {
    if (Obj==NULL) {
      if (tag==0) {
        return false;
      }
    }
    else if (!Obj->HasChanged(tag)) {
      return false;
    }
    return true;
  }

  void AugTSystemSolver::FillDPart(Index len, Number delta,
                                   const Vector* D, Number* vals2fill)
  {
    if (D) {
      // Diagonal elements are given in D
      const DenseVector* dense_d =
        dynamic_cast<const DenseVector*> (D);
      DBG_ASSERT(dense_d);
      DBG_ASSERT(dense_d->Dim()==len);
      const Number* dvals = dense_d->Values();
      IpBlasDcopy(len, dvals, 1, vals2fill, 1);
      // Now add the multiple of the identify.
      if (delta!=0.) {
        IpBlasDaxpy(len, 1.0, &delta, 0, vals2fill, 1);
      }
    }
    else {
      IpBlasDcopy(len, &delta, 0, vals2fill, 1);
    }
  }

  void AugTSystemSolver::FillWPart(Index len, const SymMatrix* W,
                                   Number* vals2fill)
  {
    if (W) {
      // Get SymTMatrix information from SymMatrix smart pointer
      const SymTMatrix* tmatrix_w =
        dynamic_cast<const SymTMatrix*>(W);
      DBG_ASSERT(tmatrix_w);
      DBG_ASSERT(tmatrix_w->Nonzeros()==len);
      const Number* wvals = tmatrix_w->Values();
      IpBlasDcopy(len, wvals, 1, vals2fill, 1);
    }
    else {
      Number zero = 0.;
      IpBlasDcopy(len, &zero, 0, vals2fill, 1);
    }
  }

  void AugTSystemSolver::FillJPart(Index len, const Matrix* J,
                                   Number* vals2fill)
  {
    if (J) {
      // Get SymTMatrix information from SymMatrix smart pointer
      const GenTMatrix* tmatrix_j =
        dynamic_cast<const GenTMatrix*>(J);
      DBG_ASSERT(tmatrix_j);
      DBG_ASSERT(tmatrix_j->Nonzeros()==len);
      const Number* jvals = tmatrix_j->Values();
      IpBlasDcopy(len, jvals, 1, vals2fill, 1);
    }
    else {
      Number zero = 0.;
      IpBlasDcopy(len, &zero, 0, vals2fill, 1);
    }
  }

  TaggedObject::Tag AugTSystemSolver::NewTag(const TaggedObject* Obj)
  {
    if (Obj==NULL) {
      return 0;
    }
    else {
      return Obj->GetTag();
    }
  }

} // namespace Ipopt
