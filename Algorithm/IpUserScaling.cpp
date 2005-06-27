// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpOrigIpoptNLP.cpp 321 2005-06-20 21:53:55Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpUserScaling.hpp"

namespace Ipopt
{

  void UserScaling::DetermineScaling(const SmartPtr<const VectorSpace> x_space,
                                     const SmartPtr<const VectorSpace> c_space,
                                     const SmartPtr<const VectorSpace> d_space,
                                     const SmartPtr<const MatrixSpace> jac_c_space,
                                     const SmartPtr<const MatrixSpace> jac_d_space,
                                     const SmartPtr<const SymMatrixSpace> h_space,
                                     SmartPtr<const MatrixSpace>& new_jac_c_space,
                                     SmartPtr<const MatrixSpace>& new_jac_d_space,
                                     SmartPtr<const SymMatrixSpace>& new_h_space)
  {
    DBG_ASSERT(IsValid(nlp_));

    dx_ = x_space->MakeNew();
    SmartPtr<Vector> dc = c_space->MakeNew();
    SmartPtr<Vector> dd = d_space->MakeNew();
    nlp_->GetScalingParameters(df_, *dx_, *dc, *dd);

    // create the scaling matrix spaces
    scaled_jac_c_space_ = new ScaledMatrixSpace(ConstPtr(dc), false, jac_c_space, ConstPtr(dx_), true);
    scaled_jac_d_space_ = new ScaledMatrixSpace(ConstPtr(dd), false, jac_d_space, ConstPtr(dx_), true);
    scaled_h_space_ = new SymScaledMatrixSpace(ConstPtr(dx_), true, h_space);

    // set the return values
    new_jac_c_space = GetRawPtr(scaled_jac_c_space_);
    new_jac_d_space = GetRawPtr(scaled_jac_d_space_);
    new_h_space = GetRawPtr(scaled_h_space_);
  }

  Number UserScaling::apply_obj_scaling(const Number& f)
  {
    return df_*f;
  }

  Number UserScaling::unapply_obj_scaling(const Number& f)
  {
    return f/df_;
  }

  SmartPtr<Vector> UserScaling::apply_vector_scaling_x_NonConst(const SmartPtr<const Vector>& v)
  {
    SmartPtr<Vector> scaled_x = v->MakeNewCopy();
    scaled_x->ElementWiseMultiply(*dx_);
    return scaled_x;
  };

  SmartPtr<const Vector> UserScaling::apply_vector_scaling_x(const SmartPtr<const Vector>& v)
  {
    return ConstPtr(apply_vector_scaling_x_NonConst(v));
  }

  SmartPtr<Vector> UserScaling::unapply_vector_scaling_x_NonConst(const SmartPtr<const Vector>& v)
  {
    SmartPtr<Vector> unscaled_x = v->MakeNewCopy();
    unscaled_x->ElementWiseDivide(*dx_);
    return unscaled_x;
  }

  SmartPtr<const Vector> UserScaling::unapply_vector_scaling_x(const SmartPtr<const Vector>& v)
  {
    return ConstPtr(unapply_vector_scaling_x_NonConst(v));
  }

  SmartPtr<Vector> UserScaling::apply_vector_scaling_c_NonConst(const SmartPtr<const Vector>& v)
  {
    SmartPtr<Vector> scaled_c = v->MakeNewCopy();
    scaled_c->ElementWiseMultiply(*scaled_jac_c_space_->RowScaling());
    return scaled_c;
  }

  SmartPtr<const Vector> UserScaling::apply_vector_scaling_c(const SmartPtr<const Vector>& v)
  {
    return ConstPtr(apply_vector_scaling_c_NonConst(v));
  }

  SmartPtr<Vector> UserScaling::unapply_vector_scaling_c_NonConst(const SmartPtr<const Vector>& v)
  {
    SmartPtr<Vector> scaled_c = v->MakeNewCopy();
    scaled_c->ElementWiseDivide(*scaled_jac_c_space_->RowScaling());
    return scaled_c;
  }

  SmartPtr<const Vector> UserScaling::unapply_vector_scaling_c(const SmartPtr<const Vector>& v)
  {
    return ConstPtr(unapply_vector_scaling_c_NonConst(v));
  }

  SmartPtr<Vector> UserScaling::apply_vector_scaling_d_NonConst(const SmartPtr<const Vector>& v)
  {
    SmartPtr<Vector> scaled_d = v->MakeNewCopy();
    scaled_d->ElementWiseMultiply(*scaled_jac_d_space_->RowScaling());
    return scaled_d;
  }

  SmartPtr<const Vector> UserScaling::apply_vector_scaling_d(const SmartPtr<const Vector>& v)
  {
    return ConstPtr(apply_vector_scaling_d_NonConst(v));
  }

  SmartPtr<Vector> UserScaling::unapply_vector_scaling_d_NonConst(const SmartPtr<const Vector>& v)
  {
    SmartPtr<Vector> scaled_d = v->MakeNewCopy();
    scaled_d->ElementWiseDivide(*scaled_jac_d_space_->RowScaling());
    return scaled_d;
  }

  SmartPtr<const Vector> UserScaling::unapply_vector_scaling_d(const SmartPtr<const Vector>& v)
  {
    return ConstPtr(unapply_vector_scaling_d_NonConst(v));
  }


  SmartPtr<const Matrix> UserScaling::apply_jac_c_scaling(SmartPtr<const Matrix> matrix)
  {
    SmartPtr<ScaledMatrix> ret = scaled_jac_c_space_->MakeNewScaledMatrix(false);
    ret->SetUnscaledMatrix(matrix);
    return GetRawPtr(ret);
  }

  SmartPtr<const Matrix> UserScaling::apply_jac_d_scaling(SmartPtr<const Matrix> matrix)
  {
    SmartPtr<ScaledMatrix> ret = scaled_jac_d_space_->MakeNewScaledMatrix(false);
    ret->SetUnscaledMatrix(matrix);
    return GetRawPtr(ret);
  }

  SmartPtr<const SymMatrix> UserScaling::apply_hessian_scaling(SmartPtr<const SymMatrix> matrix)
  {
    SmartPtr<SymScaledMatrix> ret = scaled_h_space_->MakeNewSymScaledMatrix(false);
    ret->SetUnscaledMatrix(matrix);
    return GetRawPtr(ret);
  }

} // namespace Ipopt
