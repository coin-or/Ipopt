// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpNLPScaling.hpp"

namespace Ipopt
{

  SmartPtr<Vector> NLPScalingObject::apply_vector_scaling_x_L_NonConst(SmartPtr<Matrix> Px_L, const SmartPtr<const Vector>& l, const SmartPtr<const VectorSpace> x_space)
  {
    SmartPtr<Vector> tmp_x = x_space->MakeNew();

    // move to full x space
    Px_L->MultVector(1.0, *l, 0.0, *tmp_x);

    // scale in full x space
    tmp_x = apply_vector_scaling_x_NonConst(ConstPtr(tmp_x));

    // move back to x_L space
    SmartPtr<Vector> scaled_x_L = l->MakeNew();
    Px_L->TransMultVector(1.0, *tmp_x, 0.0, *scaled_x_L);

    return scaled_x_L;
  }

  SmartPtr<Vector> NLPScalingObject::apply_vector_scaling_x_U_NonConst(SmartPtr<Matrix> Px_U, const SmartPtr<const Vector>& u, const SmartPtr<const VectorSpace> x_space)
  {
    SmartPtr<Vector> tmp_x = x_space->MakeNew();

    // move to full x space
    Px_U->MultVector(1.0, *u, 0.0, *tmp_x);

    // scale in full x space
    tmp_x = apply_vector_scaling_x_NonConst(ConstPtr(tmp_x));

    // move back to x_L space
    SmartPtr<Vector> scaled_x_U = u->MakeNew();
    Px_U->TransMultVector(1.0, *tmp_x, 0.0, *scaled_x_U);

    return scaled_x_U;
  }

  SmartPtr<Vector> NLPScalingObject::apply_vector_scaling_d_L_NonConst(SmartPtr<Matrix> Pd_L, const SmartPtr<const Vector>& l, const SmartPtr<const VectorSpace> d_space)
  {
    SmartPtr<Vector> tmp_d = d_space->MakeNew();

    // move to full d space
    Pd_L->MultVector(1.0, *l, 0.0, *tmp_d);

    // scale in full d space
    tmp_d = apply_vector_scaling_d_NonConst(ConstPtr(tmp_d));

    // move back to d_L space
    SmartPtr<Vector> scaled_d_L = l->MakeNew();
    Pd_L->TransMultVector(1.0, *tmp_d, 0.0, *scaled_d_L);

    return scaled_d_L;
  }

  SmartPtr<Vector> NLPScalingObject::apply_vector_scaling_d_U_NonConst(SmartPtr<Matrix> Pd_U, const SmartPtr<const Vector>& u, const SmartPtr<const VectorSpace> d_space)
  {
    SmartPtr<Vector> tmp_d = d_space->MakeNew();

    // move to full d space
    Pd_U->MultVector(1.0, *u, 0.0, *tmp_d);

    // scale in full d space
    tmp_d = apply_vector_scaling_d_NonConst(ConstPtr(tmp_d));

    // move back to d_L space
    SmartPtr<Vector> scaled_d_U = u->MakeNew();
    Pd_U->TransMultVector(1.0, *tmp_d, 0.0, *scaled_d_U);

    return scaled_d_U;
  }

  SmartPtr<Vector> NLPScalingObject::apply_grad_obj_scaling_NonConst(const SmartPtr<const Vector>& v)
  {
    SmartPtr<Vector> scaled_v = unapply_vector_scaling_x_NonConst(v);
    Number df = apply_obj_scaling(1.0);
    scaled_v->Scal(df);
    return scaled_v;
  }

  SmartPtr<const Vector> NLPScalingObject::apply_grad_obj_scaling(const SmartPtr<const Vector>& v)
  {
    SmartPtr<Vector> scaled_v = apply_grad_obj_scaling_NonConst(v);
    return ConstPtr(scaled_v);
  }

  SmartPtr<Vector> NLPScalingObject::unapply_grad_obj_scaling_NonConst(const SmartPtr<const Vector>& v)
  {
    SmartPtr<Vector> unscaled_v = apply_vector_scaling_x_NonConst(v);
    Number df = unapply_obj_scaling(1.0);
    unscaled_v->Scal(df);
    return unscaled_v;
  }

  SmartPtr<const Vector> NLPScalingObject::unapply_grad_obj_scaling(const SmartPtr<const Vector>& v)
  {
    SmartPtr<Vector> unscaled_v = unapply_grad_obj_scaling_NonConst(v);
    return ConstPtr(unscaled_v);
  }

  void NoNLPScalingObject::DetermineScaling(const SmartPtr<const VectorSpace> x_space,
      const SmartPtr<const VectorSpace> c_space,
      const SmartPtr<const VectorSpace> d_space,
      const SmartPtr<const MatrixSpace> jac_c_space,
      const SmartPtr<const MatrixSpace> jac_d_space,
      const SmartPtr<const SymMatrixSpace> h_space,
      SmartPtr<const MatrixSpace>& new_jac_c_space,
      SmartPtr<const MatrixSpace>& new_jac_d_space,
      SmartPtr<const SymMatrixSpace>& new_h_space)
  {
    new_jac_c_space = jac_c_space;
    new_jac_d_space = jac_d_space;
    new_h_space = h_space;
  }

  void StandardScalingBase::DetermineScaling(const SmartPtr<const VectorSpace> x_space,
      const SmartPtr<const VectorSpace> c_space,
      const SmartPtr<const VectorSpace> d_space,
      const SmartPtr<const MatrixSpace> jac_c_space,
      const SmartPtr<const MatrixSpace> jac_d_space,
      const SmartPtr<const SymMatrixSpace> h_space,
      SmartPtr<const MatrixSpace>& new_jac_c_space,
      SmartPtr<const MatrixSpace>& new_jac_d_space,
      SmartPtr<const SymMatrixSpace>& new_h_space)
  {
    dx_ = x_space->MakeNew();
    SmartPtr<Vector> dc = c_space->MakeNew();
    SmartPtr<Vector> dd = d_space->MakeNew();
    DetermineScalingParametersImpl(x_space, c_space, d_space,
                                   jac_c_space, jac_d_space,
                                   h_space,
                                   df_, *dx_, *dc, *dd);

    // create the scaling matrix spaces
    scaled_jac_c_space_ = new ScaledMatrixSpace(ConstPtr(dc), false, jac_c_space, ConstPtr(dx_), true);
    scaled_jac_d_space_ = new ScaledMatrixSpace(ConstPtr(dd), false, jac_d_space, ConstPtr(dx_), true);
    scaled_h_space_ = new SymScaledMatrixSpace(ConstPtr(dx_), true, h_space);

    // set the return values
    new_jac_c_space = GetRawPtr(scaled_jac_c_space_);
    new_jac_d_space = GetRawPtr(scaled_jac_d_space_);
    new_h_space = GetRawPtr(scaled_h_space_);
  }

  Number StandardScalingBase::apply_obj_scaling(const Number& f)
  {
    return df_*f;
  }

  Number StandardScalingBase::unapply_obj_scaling(const Number& f)
  {
    return f/df_;
  }

  SmartPtr<Vector> StandardScalingBase::apply_vector_scaling_x_NonConst(const SmartPtr<const Vector>& v)
  {
    SmartPtr<Vector> scaled_x = v->MakeNewCopy();
    scaled_x->ElementWiseMultiply(*dx_);
    return scaled_x;
  };

  SmartPtr<const Vector> StandardScalingBase::apply_vector_scaling_x(const SmartPtr<const Vector>& v)
  {
    return ConstPtr(apply_vector_scaling_x_NonConst(v));
  }

  SmartPtr<Vector> StandardScalingBase::unapply_vector_scaling_x_NonConst(const SmartPtr<const Vector>& v)
  {
    SmartPtr<Vector> unscaled_x = v->MakeNewCopy();
    unscaled_x->ElementWiseDivide(*dx_);
    return unscaled_x;
  }

  SmartPtr<const Vector> StandardScalingBase::unapply_vector_scaling_x(const SmartPtr<const Vector>& v)
  {
    return ConstPtr(unapply_vector_scaling_x_NonConst(v));
  }

  SmartPtr<Vector> StandardScalingBase::apply_vector_scaling_c_NonConst(const SmartPtr<const Vector>& v)
  {
    SmartPtr<Vector> scaled_c = v->MakeNewCopy();
    scaled_c->ElementWiseMultiply(*scaled_jac_c_space_->RowScaling());
    return scaled_c;
  }

  SmartPtr<const Vector> StandardScalingBase::apply_vector_scaling_c(const SmartPtr<const Vector>& v)
  {
    return ConstPtr(apply_vector_scaling_c_NonConst(v));
  }

  SmartPtr<Vector> StandardScalingBase::unapply_vector_scaling_c_NonConst(const SmartPtr<const Vector>& v)
  {
    SmartPtr<Vector> scaled_c = v->MakeNewCopy();
    scaled_c->ElementWiseDivide(*scaled_jac_c_space_->RowScaling());
    return scaled_c;
  }

  SmartPtr<const Vector> StandardScalingBase::unapply_vector_scaling_c(const SmartPtr<const Vector>& v)
  {
    return ConstPtr(unapply_vector_scaling_c_NonConst(v));
  }

  SmartPtr<Vector> StandardScalingBase::apply_vector_scaling_d_NonConst(const SmartPtr<const Vector>& v)
  {
    SmartPtr<Vector> scaled_d = v->MakeNewCopy();
    scaled_d->ElementWiseMultiply(*scaled_jac_d_space_->RowScaling());
    return scaled_d;
  }

  SmartPtr<const Vector> StandardScalingBase::apply_vector_scaling_d(const SmartPtr<const Vector>& v)
  {
    return ConstPtr(apply_vector_scaling_d_NonConst(v));
  }

  SmartPtr<Vector> StandardScalingBase::unapply_vector_scaling_d_NonConst(const SmartPtr<const Vector>& v)
  {
    SmartPtr<Vector> scaled_d = v->MakeNewCopy();
    scaled_d->ElementWiseDivide(*scaled_jac_d_space_->RowScaling());
    return scaled_d;
  }

  SmartPtr<const Vector> StandardScalingBase::unapply_vector_scaling_d(const SmartPtr<const Vector>& v)
  {
    return ConstPtr(unapply_vector_scaling_d_NonConst(v));
  }


  SmartPtr<const Matrix> StandardScalingBase::apply_jac_c_scaling(SmartPtr<const Matrix> matrix)
  {
    SmartPtr<ScaledMatrix> ret = scaled_jac_c_space_->MakeNewScaledMatrix(false);
    ret->SetUnscaledMatrix(matrix);
    return GetRawPtr(ret);
  }

  SmartPtr<const Matrix> StandardScalingBase::apply_jac_d_scaling(SmartPtr<const Matrix> matrix)
  {
    SmartPtr<ScaledMatrix> ret = scaled_jac_d_space_->MakeNewScaledMatrix(false);
    ret->SetUnscaledMatrix(matrix);
    return GetRawPtr(ret);
  }

  SmartPtr<const SymMatrix> StandardScalingBase::apply_hessian_scaling(SmartPtr<const SymMatrix> matrix)
  {
    SmartPtr<SymScaledMatrix> ret = scaled_h_space_->MakeNewSymScaledMatrix(false);
    ret->SetUnscaledMatrix(matrix);
    return GetRawPtr(ret);
  }

} // namespace Ipopt
