// Copyright (C) 2008, 2010 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: IpTNLP.hpp 1235 2008-05-22 14:38:40Z andreasw $
//
// Authors:  Andreas Waechter                  IBM    2008-08-25

#include "IpNLPBoundsRemover.hpp"
#include "IpCompoundVector.hpp"
#include "IpCompoundMatrix.hpp"
#include "IpCompoundSymMatrix.hpp"
#include "IpIdentityMatrix.hpp"
#include "IpTransposeMatrix.hpp"
#include "IpDenseVector.hpp"
#include "IpZeroMatrix.hpp"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  NLPBoundsRemover::NLPBoundsRemover(NLP& nlp,
                                     bool allow_twosided_inequalities /* = false */)
      :
      nlp_(&nlp),
      allow_twosided_inequalities_(allow_twosided_inequalities)
  {}

  bool
  NLPBoundsRemover::GetSpaces(SmartPtr<const VectorSpace>& x_space,
                              SmartPtr<const VectorSpace>& c_space,
                              SmartPtr<const VectorSpace>& d_space,
                              SmartPtr<const VectorSpace>& x_l_space,
                              SmartPtr<const MatrixSpace>& px_l_space,
                              SmartPtr<const VectorSpace>& x_u_space,
                              SmartPtr<const MatrixSpace>& px_u_space,
                              SmartPtr<const VectorSpace>& d_l_space,
                              SmartPtr<const MatrixSpace>& pd_l_space,
                              SmartPtr<const VectorSpace>& d_u_space,
                              SmartPtr<const MatrixSpace>& pd_u_space,
                              SmartPtr<const MatrixSpace>& Jac_c_space,
                              SmartPtr<const MatrixSpace>& Jac_d_space,
                              SmartPtr<const SymMatrixSpace>& Hess_lagrangian_space)
  {
    DBG_START_METH("NLPBoundsRemover::GetSpaces", dbg_verbosity);
    SmartPtr<const VectorSpace> d_space_orig;
    SmartPtr<const VectorSpace> x_l_space_orig;
    SmartPtr<const MatrixSpace> px_l_space_orig;
    SmartPtr<const VectorSpace> x_u_space_orig;
    SmartPtr<const MatrixSpace> px_u_space_orig;
    SmartPtr<const VectorSpace> d_l_space_orig;
    SmartPtr<const MatrixSpace> pd_l_space_orig;
    SmartPtr<const VectorSpace> d_u_space_orig;
    SmartPtr<const MatrixSpace> pd_u_space_orig;
    SmartPtr<const MatrixSpace> Jac_d_space_orig;

    bool retval = nlp_->GetSpaces(x_space, c_space, d_space_orig,
                                  x_l_space_orig, px_l_space_orig,
                                  x_u_space_orig, px_u_space_orig,
                                  d_l_space_orig, pd_l_space_orig,
                                  d_u_space_orig, pd_u_space_orig,
                                  Jac_c_space, Jac_d_space_orig,
                                  Hess_lagrangian_space);
    if (!retval) {
      return retval;
    }
    // Keep a copy of the expansion matrices for the x bounds
    Px_l_orig_ = px_l_space_orig->MakeNew();
    Px_u_orig_ = px_u_space_orig->MakeNew();

    // create the new d_space
    Index total_dim = d_space_orig->Dim() + x_l_space_orig->Dim() +
                      x_u_space_orig->Dim();
    SmartPtr<CompoundVectorSpace> d_space_new =
      new CompoundVectorSpace(3, total_dim);
    d_space_new->SetCompSpace(0, *d_space_orig);
    d_space_new->SetCompSpace(1, *x_l_space_orig);
    d_space_new->SetCompSpace(2, *x_u_space_orig);
    d_space = GetRawPtr(d_space_new);

    // create the new (emply) x_l and x_u spaces, and also the
    // corresponding projection matrix spaces
    x_l_space = new DenseVectorSpace(0);
    x_u_space = new DenseVectorSpace(0);
    px_l_space = new ZeroMatrixSpace(x_space->Dim(),0);
    px_u_space = new ZeroMatrixSpace(x_space->Dim(),0);

    // create the new d_l and d_u vector spaces
    total_dim = d_l_space_orig->Dim() + x_l_space_orig->Dim();
    SmartPtr<CompoundVectorSpace> d_l_space_new =
      new CompoundVectorSpace(2, total_dim);
    d_l_space_new->SetCompSpace(0, *d_l_space_orig);
    d_l_space_new->SetCompSpace(1, *x_l_space_orig);
    d_l_space = GetRawPtr(d_l_space_new);
    total_dim = d_u_space_orig->Dim() + x_u_space_orig->Dim();
    SmartPtr<CompoundVectorSpace> d_u_space_new =
      new CompoundVectorSpace(2, total_dim);
    d_u_space_new->SetCompSpace(0, *d_u_space_orig);
    d_u_space_new->SetCompSpace(1, *x_u_space_orig);
    d_u_space = GetRawPtr(d_u_space_new);

    // create the new d_l and d_u projection matrix spaces
    Index total_rows = d_space_orig->Dim() + x_l_space_orig->Dim() +
                       x_u_space_orig->Dim();
    Index total_cols = d_l_space_orig->Dim() + x_l_space_orig->Dim();
    SmartPtr<CompoundMatrixSpace> pd_l_space_new =
      new CompoundMatrixSpace(3, 2, total_rows, total_cols);
    pd_l_space_new->SetBlockRows(0, d_space_orig->Dim());
    pd_l_space_new->SetBlockRows(1, x_l_space_orig->Dim());
    pd_l_space_new->SetBlockRows(2, x_u_space_orig->Dim());
    pd_l_space_new->SetBlockCols(0, d_l_space_orig->Dim());
    pd_l_space_new->SetBlockCols(1, x_l_space_orig->Dim());
    pd_l_space_new->SetCompSpace(0, 0, *pd_l_space_orig, true);
    SmartPtr<const MatrixSpace> identity_space =
      new IdentityMatrixSpace(x_l_space_orig->Dim());
    pd_l_space_new->SetCompSpace(1, 1, *identity_space, true);
    pd_l_space = GetRawPtr(pd_l_space_new);

    total_cols = d_u_space_orig->Dim() + x_u_space_orig->Dim();
    SmartPtr<CompoundMatrixSpace> pd_u_space_new =
      new CompoundMatrixSpace(3, 2, total_rows, total_cols);
    pd_u_space_new->SetBlockRows(0, d_space_orig->Dim());
    pd_u_space_new->SetBlockRows(1, x_l_space_orig->Dim());
    pd_u_space_new->SetBlockRows(2, x_u_space_orig->Dim());
    pd_u_space_new->SetBlockCols(0, d_u_space_orig->Dim());
    pd_u_space_new->SetBlockCols(1, x_u_space_orig->Dim());
    pd_u_space_new->SetCompSpace(0, 0, *pd_u_space_orig, true);
    identity_space = new IdentityMatrixSpace(x_u_space_orig->Dim());
    pd_u_space_new->SetCompSpace(2, 1, *identity_space, true);
    pd_u_space = GetRawPtr(pd_u_space_new);

    // Jacobian for inequalities matrix space
    total_rows = d_space_orig->Dim() + x_l_space_orig->Dim() +
                 x_u_space_orig->Dim();
    total_cols = x_space->Dim();
    SmartPtr<CompoundMatrixSpace> Jac_d_space_new =
      new CompoundMatrixSpace(3, 1, total_rows, total_cols);
    Jac_d_space_new->SetBlockRows(0, d_space_orig->Dim());
    Jac_d_space_new->SetBlockRows(1, x_l_space_orig->Dim());
    Jac_d_space_new->SetBlockRows(2, x_u_space_orig->Dim());
    Jac_d_space_new->SetBlockCols(0, x_space->Dim());
    Jac_d_space_new->SetCompSpace(0, 0, *Jac_d_space_orig);
    SmartPtr<MatrixSpace> trans_px_l_space_orig =
      new TransposeMatrixSpace(GetRawPtr(px_l_space_orig));
    Jac_d_space_new->SetCompSpace(1, 0, *trans_px_l_space_orig, true);
    SmartPtr<MatrixSpace> trans_px_u_space_orig =
      new TransposeMatrixSpace(GetRawPtr(px_u_space_orig));
    Jac_d_space_new->SetCompSpace(2, 0, *trans_px_u_space_orig, true);
    Jac_d_space = GetRawPtr(Jac_d_space_new);

    // We keep the original d_space around in order to be able to do
    // the sanity check later
    d_space_orig_ = d_space_orig;

    return true;
  }

  bool
  NLPBoundsRemover::GetBoundsInformation(const Matrix& Px_L,
                                         Vector& x_L,
                                         const Matrix& Px_U,
                                         Vector& x_U,
                                         const Matrix& Pd_L,
                                         Vector& d_L,
                                         const Matrix& Pd_U,
                                         Vector& d_U)
  {
    const CompoundMatrix* comp_pd_l =
      static_cast<const CompoundMatrix*>(&Pd_L);
    DBG_ASSERT(dynamic_cast<const CompoundMatrix*>(&Pd_L));
    SmartPtr<const Matrix> pd_l_orig = comp_pd_l->GetComp(0,0);

    const CompoundMatrix* comp_pd_u =
      static_cast<const CompoundMatrix*>(&Pd_U);
    DBG_ASSERT(dynamic_cast<const CompoundMatrix*>(&Pd_U));
    SmartPtr<const Matrix> pd_u_orig = comp_pd_u->GetComp(0,0);

    CompoundVector* comp_d_l = static_cast<CompoundVector*>(&d_L);
    DBG_ASSERT(dynamic_cast<CompoundVector*>(&d_L));
    SmartPtr<Vector> d_l_orig = comp_d_l->GetCompNonConst(0);
    SmartPtr<Vector> x_l_orig = comp_d_l->GetCompNonConst(1);

    CompoundVector* comp_d_u = static_cast<CompoundVector*>(&d_U);
    DBG_ASSERT(dynamic_cast<CompoundVector*>(&d_U));
    SmartPtr<Vector> d_u_orig = comp_d_u->GetCompNonConst(0);
    SmartPtr<Vector> x_u_orig = comp_d_u->GetCompNonConst(1);

    // Here we do a santiy check to make sure that no inequality
    // constraint has two non-infite bounds.
    if (d_space_orig_->Dim()>0 && !allow_twosided_inequalities_) {
      SmartPtr<Vector> d = d_space_orig_->MakeNew();
      SmartPtr<Vector> tmp = d_l_orig->MakeNew();
      tmp->Set(1.);
      pd_l_orig->MultVector(1., *tmp, 0., *d);
      tmp = d_u_orig->MakeNew();
      tmp->Set(1.);
      pd_u_orig->MultVector(1., *tmp, 1., *d);
      Number dmax = d->Amax();
      ASSERT_EXCEPTION(dmax==1., INVALID_NLP, "In NLPBoundRemover, an inequality with both lower and upper bounds was detected");
      Number dmin = d->Min();
      ASSERT_EXCEPTION(dmin==1., INVALID_NLP, "In NLPBoundRemover, an inequality with without bounds was detected.");
    }

    bool retval =
      nlp_->GetBoundsInformation(*Px_l_orig_, *x_l_orig, *Px_u_orig_,
                                 *x_u_orig, *pd_l_orig, *d_l_orig,
                                 *pd_u_orig, *d_u_orig);
    return retval;
  }

  bool
  NLPBoundsRemover::GetStartingPoint(SmartPtr<Vector> x,
                                     bool need_x,
                                     SmartPtr<Vector> y_c,
                                     bool need_y_c,
                                     SmartPtr<Vector> y_d,
                                     bool need_y_d,
                                     SmartPtr<Vector> z_L,
                                     bool need_z_L,
                                     SmartPtr<Vector> z_U,
                                     bool need_z_U)
  {
    SmartPtr<Vector> y_d_orig;
    SmartPtr<Vector> z_L_orig;
    SmartPtr<Vector> z_U_orig;
    if (need_y_d) {
      CompoundVector* comp_y_d = static_cast<CompoundVector*>(GetRawPtr(y_d));
      DBG_ASSERT(dynamic_cast<CompoundVector*>(GetRawPtr(y_d)));
      y_d_orig = comp_y_d->GetCompNonConst(0);
      z_L_orig = comp_y_d->GetCompNonConst(1);
      z_U_orig = comp_y_d->GetCompNonConst(2);
    }
    bool retval =
      nlp_->GetStartingPoint(x, need_x, y_c, need_y_c, y_d_orig, need_y_d,
                             z_L_orig, need_y_d, z_U_orig, need_y_d);
    return retval;
  }

  bool
  NLPBoundsRemover::Eval_d(const Vector& x, Vector& d)
  {
    CompoundVector* comp_d = static_cast<CompoundVector*>(&d);
    DBG_ASSERT(dynamic_cast<CompoundVector*>(&d));
    SmartPtr<Vector> d_orig = comp_d->GetCompNonConst(0);

    bool retval = nlp_->Eval_d(x, *d_orig);
    if (retval) {
      SmartPtr<Vector> x_L = comp_d->GetCompNonConst(1);
      SmartPtr<Vector> x_U = comp_d->GetCompNonConst(2);
      Px_l_orig_->TransMultVector(1., x, 0., *x_L);
      Px_u_orig_->TransMultVector(1., x, 0., *x_U);
    }
    return retval;
  }

  bool
  NLPBoundsRemover::Eval_jac_d(const Vector& x, Matrix& jac_d)
  {
    CompoundMatrix* comp_jac_d = static_cast<CompoundMatrix*>(&jac_d);
    DBG_ASSERT(dynamic_cast<CompoundMatrix*>(&jac_d));
    SmartPtr<const MatrixSpace> jac_d_space = comp_jac_d->OwnerSpace();
    const CompoundMatrixSpace* comp_jac_d_space =
      static_cast<const CompoundMatrixSpace*>(GetRawPtr(jac_d_space));
    DBG_ASSERT(dynamic_cast<const CompoundMatrixSpace*>(GetRawPtr(jac_d_space)));
    SmartPtr<Matrix> jac_d_orig = comp_jac_d_space->GetCompSpace(0,0)->MakeNew();
    bool retval = nlp_->Eval_jac_d(x, *jac_d_orig);
    if (retval) {
      comp_jac_d->SetComp(0, 0, *jac_d_orig);
    }
    return retval;
  }

  bool
  NLPBoundsRemover::Eval_h(const Vector& x, Number obj_factor,
                           const Vector& yc, const Vector& yd, SymMatrix& h)
  {
    const CompoundVector* comp_yd = static_cast<const CompoundVector*>(&yd);
    DBG_ASSERT(dynamic_cast<const CompoundVector*>(&yd));
    SmartPtr<const Vector> yd_orig = comp_yd->GetComp(0);

    bool retval = nlp_->Eval_h(x, obj_factor, yc, *yd_orig, h);
    return retval;
  }

  void
  NLPBoundsRemover::FinalizeSolution(SolverReturn status,
                                     const Vector& x, const Vector& z_L,
                                     const Vector& z_U,
                                     const Vector& c, const Vector& d,
                                     const Vector& y_c, const Vector& y_d,
                                     Number obj_value,
                                     const IpoptData* ip_data,
                                     IpoptCalculatedQuantities* ip_cq)
  {
    const CompoundVector* comp_d = static_cast<const CompoundVector*>(&d);
    DBG_ASSERT(dynamic_cast<const CompoundVector*>(&d));
    SmartPtr<const Vector> d_orig = comp_d->GetComp(0);

    const CompoundVector* comp_y_d = static_cast<const CompoundVector*>(&y_d);
    DBG_ASSERT(dynamic_cast<const CompoundVector*>(&y_d));
    SmartPtr<const Vector> y_d_orig = comp_y_d->GetComp(0);
    SmartPtr<const Vector> z_L_orig = comp_y_d->GetComp(1);
    SmartPtr<const Vector> z_U_orig = comp_y_d->GetComp(2);

    SmartPtr<Vector> z_L_new = z_L_orig->MakeNewCopy();
    z_L_new->Scal(-1.);

    nlp_->FinalizeSolution(status, x, *z_L_new, *z_U_orig, c, *d_orig,
                           y_c, *y_d_orig, obj_value, ip_data, ip_cq);
  }

  void
  NLPBoundsRemover::GetScalingParameters(
    const SmartPtr<const VectorSpace> x_space,
    const SmartPtr<const VectorSpace> c_space,
    const SmartPtr<const VectorSpace> d_space,
    Number& obj_scaling,
    SmartPtr<Vector>& x_scaling,
    SmartPtr<Vector>& c_scaling,
    SmartPtr<Vector>& d_scaling) const
  {
    const CompoundVectorSpace* comp_d_space =
      static_cast<const CompoundVectorSpace*>(GetRawPtr(d_space));
    DBG_ASSERT(dynamic_cast<const CompoundVectorSpace*>(GetRawPtr(d_space)));
    SmartPtr<const VectorSpace> d_space_orig = comp_d_space->GetCompSpace(0);

    SmartPtr<Vector> d_scaling_orig;
    nlp_->GetScalingParameters(x_space, c_space, d_space_orig, obj_scaling,
                               x_scaling, c_scaling, d_scaling_orig);

    if (IsValid(x_scaling) || IsValid(d_scaling_orig)) {

      SmartPtr<CompoundVector> comp_d_scaling =
        comp_d_space->MakeNewCompoundVector();

      SmartPtr<Vector> xL_scaling = comp_d_scaling->GetCompNonConst(1);
      SmartPtr<Vector> xU_scaling = comp_d_scaling->GetCompNonConst(2);
      if (IsValid(x_scaling)) {
        Px_l_orig_->TransMultVector(1., *x_scaling, 0., *xL_scaling);
        Px_u_orig_->TransMultVector(1., *x_scaling, 0., *xU_scaling);
      }
      else {
        xL_scaling->Set(1.);
        xU_scaling->Set(1.);
      }

      if (IsValid(d_scaling_orig)) {
        comp_d_scaling->SetComp(0, *d_scaling_orig);
      }
      else {
        comp_d_scaling->GetCompNonConst(0)->Set(1.);
      }

      d_scaling = GetRawPtr(comp_d_scaling);
    }
    else {
      d_scaling = NULL;
    }
  }

}

