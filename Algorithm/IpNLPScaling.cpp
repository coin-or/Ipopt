// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpOrigIpoptNLP.cpp 321 2005-06-20 21:53:55Z andreasw $
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
} // namespace Ipopt
