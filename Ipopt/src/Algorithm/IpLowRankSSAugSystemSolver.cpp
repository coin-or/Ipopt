// Copyright (C) 2005, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter                IBM    2009-11-05
//             (based on IpLowRankAugSystemSolver.cpp rev 1571)

#include "IpLowRankSSAugSystemSolver.hpp"
#include "IpLowRankUpdateSymMatrix.hpp"
#include "IpDenseVector.hpp"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  LowRankSSAugSystemSolver::LowRankSSAugSystemSolver(
    AugSystemSolver& aug_system_solver,
    Index max_rank)
      :
      AugSystemSolver(),
      aug_system_solver_(&aug_system_solver),
      max_rank_(max_rank),
      w_tag_(0),
      w_factor_(0.),
      d_x_tag_(0),
      delta_x_(0.),
      d_s_tag_(0),
      delta_s_(0.),
      j_c_tag_(0),
      d_c_tag_(0),
      delta_c_(0.),
      j_d_tag_(0),
      d_d_tag_(0),
      delta_d_(0.)
  {
    DBG_START_METH("LowRankSSAugSystemSolver::LowRankSSAugSystemSolver()",dbg_verbosity);
    DBG_ASSERT(IsValid(aug_system_solver_));
  }

  LowRankSSAugSystemSolver::~LowRankSSAugSystemSolver()
  {
    DBG_START_METH("LowRankSSAugSystemSolver::~LowRankSSAugSystemSolver()",dbg_verbosity);
  }

  bool LowRankSSAugSystemSolver::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    first_call_ = true;
    Wdiag_ = NULL;
    expanded_vu_ = NULL;
    J_c_ext_ = NULL;
    D_c_ext_ = NULL;
    y_c_ext_space_ = NULL;

    return aug_system_solver_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                          options, prefix);
  }

  ESymSolverStatus LowRankSSAugSystemSolver::Solve(
    const SymMatrix* W,
    double W_factor,
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
    DBG_START_METH("LowRankSSAugSystemSolver::Solve",dbg_verbosity);

    ESymSolverStatus retval;

    if (first_call_) {
      DBG_ASSERT(IsNull(Wdiag_));
      // Set up the diagonal matrix Wdiag_
      Index dimx = rhs_x.Dim();
      SmartPtr<DiagMatrixSpace> Wdiag_space = new DiagMatrixSpace(dimx);
      Wdiag_ = Wdiag_space->MakeNewDiagMatrix();
    }

    // This might be used with a linear solver that cannot detect the
    // inertia.  In that case, we should not asked for checking the
    // number of negative eigenvalues.
    if (!aug_system_solver_->ProvidesInertia()) {
      check_NegEVals = false;
    }

    if (first_call_ ||
        AugmentedSystemRequiresChange(W, W_factor, D_x, delta_x, D_s, delta_s,
                                      *J_c, D_c, delta_c, *J_d, D_d,
                                      delta_d) ) {
      retval = UpdateExtendedData(W, W_factor, D_x, delta_x, D_s, delta_s,
                                  *J_c, D_c, delta_c, *J_d, D_d, delta_d,
                                  rhs_x, rhs_s, rhs_c, rhs_d);
      if (retval != SYMSOLVER_SUCCESS) {
        return retval;
      }

      // Store the tags
      w_tag_ = W->GetTag();
      w_factor_ = W_factor;
      if (D_x) {
        d_x_tag_ = D_x->GetTag();
      }
      else {
        d_x_tag_ = 0;
      }
      delta_x_ = delta_x;
      if (D_s) {
        d_s_tag_ = D_s->GetTag();
      }
      else {
        d_s_tag_ = 0;
      }
      delta_s_ = delta_s;
      if (J_c) {
        j_c_tag_ = J_c->GetTag();
      }
      else {
        j_c_tag_ = 0;
      }
      if (D_c) {
        d_c_tag_ = D_c->GetTag();
      }
      else {
        d_c_tag_ = 0;
      }
      delta_c_ = delta_c;
      if (J_d) {
        j_d_tag_ = J_d->GetTag();
      }
      else {
        j_d_tag_ = 0;
      }
      if (D_d) {
        d_d_tag_ = D_d->GetTag();
      }
      else {
        d_d_tag_ = 0;
      }
      delta_d_ = delta_d;

      first_call_ = false;
    }

    // Extend the right hand side
    SmartPtr<CompoundVector> rhs_c_ext =
      y_c_ext_space_->MakeNewCompoundVector(true);
    rhs_c_ext->SetComp(0, rhs_c);
    rhs_c_ext->GetCompNonConst(1)->Set(0.);
    SmartPtr<CompoundVector> sol_c_ext =
      y_c_ext_space_->MakeNewCompoundVector(true);
    sol_c_ext->SetCompNonConst(0, sol_c);

    // Now solve the system for the given right hand side, using the
    // extended Jacobian_c and y_c data.
    numberOfNegEVals += negEvalsCorrection_;
    retval = aug_system_solver_->Solve(GetRawPtr(Wdiag_), W_factor,
                                       D_x, delta_x, D_s, delta_s,
                                       GetRawPtr(J_c_ext_),
                                       GetRawPtr(D_c_ext_), delta_c,
                                       J_d, D_d, delta_d,
                                       rhs_x, rhs_s, *rhs_c_ext, rhs_d,
                                       sol_x, sol_s, *sol_c_ext, sol_d,
                                       check_NegEVals, numberOfNegEVals);
    if (aug_system_solver_->ProvidesInertia()) {
      num_neg_evals_ =
        aug_system_solver_->NumberOfNegEVals() - negEvalsCorrection_;
    }
    if (retval != SYMSOLVER_SUCCESS) {
      Jnlst().Printf(J_DETAILED, J_SOLVE_PD_SYSTEM,
                     "LowRankSSAugSystemSolver: AugSystemSolver returned retval = %d for right hand side.\n", retval);
      return retval;
    }

    return retval;
  }

  ESymSolverStatus LowRankSSAugSystemSolver::UpdateExtendedData(
    const SymMatrix* W,
    double W_factor,
    const Vector* D_x,
    double delta_x,
    const Vector* D_s,
    double delta_s,
    const Matrix& J_c,
    const Vector* D_c,
    double delta_c,
    const Matrix& J_d,
    const Vector* D_d,
    double delta_d,
    const Vector& proto_rhs_x,
    const Vector& proto_rhs_s,
    const Vector& proto_rhs_c,
    const Vector& proto_rhs_d)
  {
    DBG_START_METH("LowRankSSAugSystemSolver::UpdateExtendedData",
                   dbg_verbosity);

    DBG_ASSERT(W_factor == 0.0 || W_factor == 1.0);
    ESymSolverStatus retval = SYMSOLVER_SUCCESS;

    // Get the low update information out of W
    const LowRankUpdateSymMatrix* LR_W =
      static_cast<const LowRankUpdateSymMatrix*> (W);
    DBG_ASSERT(dynamic_cast<const LowRankUpdateSymMatrix*>(W));
    DBG_PRINT_MATRIX(2, "LR_W", *LR_W);

    // If we don't have it yet, create the ExpandedMultiVectorMatrix
    SmartPtr<const Matrix> P_LM = LR_W->P_LowRank();
    SmartPtr<const VectorSpace> LR_VecSpace = LR_W->LowRankVectorSpace();
    if (IsNull(expanded_vu_)) {
      SmartPtr<const ExpansionMatrix> exp_matrix;
      if (IsValid(P_LM)) {
        exp_matrix = static_cast<const ExpansionMatrix*>(GetRawPtr(P_LM));
        DBG_ASSERT(dynamic_cast<const ExpansionMatrix*>(GetRawPtr(P_LM)));
      }
      SmartPtr<ExpandedMultiVectorMatrixSpace> expanded_vu_space =
        new ExpandedMultiVectorMatrixSpace(max_rank_, *LR_VecSpace, exp_matrix);
      expanded_vu_ = expanded_vu_space->MakeNewExpandedMultiVectorMatrix();

      // Create extended y_c quantities to include the V and U matrices
      DBG_ASSERT(IsNull(J_c_ext_));
      SmartPtr<CompoundMatrixSpace> J_c_ext_space =
        new CompoundMatrixSpace(2, 1, proto_rhs_c.Dim()+max_rank_,
                                proto_rhs_x.Dim());
      J_c_ext_space->SetBlockRows(0, proto_rhs_c.Dim());
      J_c_ext_space->SetBlockRows(1, max_rank_);
      J_c_ext_space->SetBlockCols(0, proto_rhs_x.Dim());
      J_c_ext_space->SetCompSpace(0, 0, *J_c.OwnerSpace());
      J_c_ext_space->SetCompSpace(1, 0, *expanded_vu_space);

      J_c_ext_ = J_c_ext_space->MakeNewCompoundMatrix();

      DBG_ASSERT(IsNull(D_c_ext_));
      DBG_ASSERT(IsNull(y_c_ext_space_));
      y_c_ext_space_ = new CompoundVectorSpace(2, proto_rhs_c.Dim()+max_rank_);
      y_c_ext_space_->SetCompSpace(0, *proto_rhs_c.OwnerSpace());
      SmartPtr<DenseVectorSpace> D_c_rank_space =
        new DenseVectorSpace(max_rank_);
      y_c_ext_space_->SetCompSpace(1, *D_c_rank_space);
      D_c_ext_ = y_c_ext_space_->MakeNewCompoundVector(true);
    }

    SmartPtr<const Vector> B0;
    SmartPtr<const MultiVectorMatrix> V;
    SmartPtr<const MultiVectorMatrix> U;
    if (W_factor == 1.0) {
      V = LR_W->GetV();
      U = LR_W->GetU();
      B0 = LR_W->GetDiag();
    }

    if (IsNull(B0)) {
      SmartPtr<Vector> zero_B0 = (IsValid(P_LM)) ? LR_VecSpace->MakeNew() : proto_rhs_x.MakeNew();
      zero_B0->Set(0.0);
      B0 = GetRawPtr(zero_B0);
    }

    // set up the Hessian for the underlying augmented system solver
    // without the low-rank update
    if (IsValid(P_LM) && LR_W->ReducedDiag()) {
      DBG_ASSERT(IsValid(B0));
      SmartPtr<Vector> fullx = proto_rhs_x.MakeNew();
      P_LM->MultVector(1., *B0, 0., *fullx);
      Wdiag_->SetDiag(*fullx);
    }
    else {
      Wdiag_->SetDiag(*B0);
      DBG_PRINT_VECTOR(2, "B0", *B0);
    }

    SmartPtr<Vector> D_c_rank_vec = D_c_ext_->GetCompNonConst(1);
    SmartPtr<DenseVector> D_c_rank =
      static_cast<DenseVector*>(GetRawPtr(D_c_rank_vec));
    DBG_ASSERT(dynamic_cast<DenseVector*>(GetRawPtr(D_c_rank_vec)));
    Number* D_c_rank_vals = D_c_rank->Values();
    Index irank = 0;
    if (IsValid(V)) {
      Index nV = V->NCols();
      negEvalsCorrection_ = nV;
      ASSERT_EXCEPTION(irank + nV, INTERNAL_ABORT, "max_rank too small for V");
      for (Index i=0; i<nV; i++) {
        SmartPtr<const Vector> vec = V->GetVector(i);
        expanded_vu_->SetVector(irank, vec);
        D_c_rank_vals[irank] = -1.;
        irank++;
      }
    }
    else {
      negEvalsCorrection_ = 0;
    }
    if (IsValid(U)) {
      Index nU = U->NCols();
      ASSERT_EXCEPTION(irank + nU, INTERNAL_ABORT, "max_rank too small for V");
      for (Index i=0; i<nU; i++) {
        SmartPtr<const Vector> vec = U->GetVector(i);
        expanded_vu_->SetVector(irank, vec);
        D_c_rank_vals[irank] = 1.;
        irank++;
      }
    }
    for (; irank<max_rank_; irank++) {
      expanded_vu_->SetVector(irank, NULL);
      D_c_rank_vals[irank] = 1.;
    }
    if (D_c) {
      D_c_ext_->SetComp(0, *D_c);
    }
    else {
      SmartPtr<Vector> zero_c = proto_rhs_c.MakeNew();
      zero_c->Set(0.);
      D_c_ext_->SetComp(0, *zero_c);
    }
    J_c_ext_->SetComp(0, 0, J_c);
    J_c_ext_->SetComp(1, 0, *expanded_vu_);

    return retval;
  }

  bool LowRankSSAugSystemSolver::AugmentedSystemRequiresChange(
    const SymMatrix* W,
    double W_factor,
    const Vector* D_x,
    double delta_x,
    const Vector* D_s,
    double delta_s,
    const Matrix& J_c,
    const Vector* D_c,
    double delta_c,
    const Matrix& J_d,
    const Vector* D_d,
    double delta_d)
  {
    DBG_START_METH("LowRankSSAugSystemSolver::AugmentedSystemRequiresChange",
                   dbg_verbosity);

#if COIN_IPOPT_VERBOSITY > 0

    bool Wtest = (W && W->GetTag() != w_tag_);
    bool iWtest = (!W && w_tag_ != 0);
    bool wfactor_test = (W_factor != w_factor_);
    bool D_xtest = (D_x && D_x->GetTag() != d_x_tag_);
    bool iD_xtest = (!D_x && d_x_tag_ != 0);
    bool delta_xtest = (delta_x != delta_x_);
    bool D_stest = (D_s && D_s->GetTag() != d_s_tag_);
    bool iD_stest = (!D_s && d_s_tag_ != 0);
    bool delta_stest = (delta_s != delta_s_);
    bool J_ctest = (J_c.GetTag() != j_c_tag_);
    bool D_ctest = (D_c && D_c->GetTag() != d_c_tag_);
    bool iD_ctest = (!D_c && d_c_tag_ != 0);
    bool delta_ctest = (delta_c != delta_c_);
    bool J_dtest = (J_d.GetTag() != j_d_tag_);
    bool D_dtest = (D_d && D_d->GetTag() != d_d_tag_);
    bool iD_dtest = (!D_d && d_d_tag_ != 0);
    bool delta_dtest = (delta_d != delta_d_);
#endif

    DBG_PRINT((2,"Wtest = %d\n", Wtest));
    DBG_PRINT((2,"iWtest = %d\n", iWtest));
    DBG_PRINT((2,"wfactor_test = %d\n", wfactor_test));
    DBG_PRINT((2,"D_xtest = %d\n", D_xtest));
    DBG_PRINT((2,"iD_xtest = %d\n", iD_xtest));
    DBG_PRINT((2,"delta_xtest = %d\n", delta_xtest));
    DBG_PRINT((2,"D_stest = %d\n", D_stest));
    DBG_PRINT((2,"iD_stest = %d\n", iD_stest));
    DBG_PRINT((2,"delta_stest = %d\n", delta_stest));
    DBG_PRINT((2,"J_ctest = %d\n", J_ctest));
    DBG_PRINT((2,"D_ctest = %d\n", D_ctest));
    DBG_PRINT((2,"iD_ctest = %d\n", iD_ctest));
    DBG_PRINT((2,"delta_ctest = %d\n", delta_ctest));
    DBG_PRINT((2,"J_dtest = %d\n", J_dtest));
    DBG_PRINT((2,"D_dtest = %d\n", D_dtest));
    DBG_PRINT((2,"iD_dtest = %d\n", iD_dtest));
    DBG_PRINT((2,"delta_dtest = %d\n", delta_dtest));

    if ( (W && W->GetTag() != w_tag_)
         || (!W && w_tag_ != 0)
         || (W_factor != w_factor_)
         || (D_x && D_x->GetTag() != d_x_tag_)
         || (!D_x && d_x_tag_ != 0)
         || (delta_x != delta_x_)
         || (D_s && D_s->GetTag() != d_s_tag_)
         || (!D_s && d_s_tag_ != 0)
         || (delta_s != delta_s_)
         || (J_c.GetTag() != j_c_tag_)
         || (D_c && D_c->GetTag() != d_c_tag_)
         || (!D_c && d_c_tag_ != 0)
         || (delta_c != delta_c_)
         || (J_d.GetTag() != j_d_tag_)
         || (D_d && D_d->GetTag() != d_d_tag_)
         || (!D_d && d_d_tag_ != 0)
         || (delta_d != delta_d_) ) {
      return true;
    }

    return false;
  }

  Index LowRankSSAugSystemSolver::NumberOfNegEVals() const
  {
    DBG_ASSERT(!first_call_);
    return num_neg_evals_;
  }

  bool LowRankSSAugSystemSolver::ProvidesInertia() const
  {
    return aug_system_solver_->ProvidesInertia();
  }

  bool LowRankSSAugSystemSolver::IncreaseQuality()
  {
    return aug_system_solver_->IncreaseQuality();
  }

} // namespace Ipopt
