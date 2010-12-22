// Copyright (C) 2004, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-06-25

#include "IpUserScaling.hpp"

namespace Ipopt
{

  void UserScaling::DetermineScalingParametersImpl(
    const SmartPtr<const VectorSpace> x_space,
    const SmartPtr<const VectorSpace> c_space,
    const SmartPtr<const VectorSpace> d_space,
    const SmartPtr<const MatrixSpace> jac_c_space,
    const SmartPtr<const MatrixSpace> jac_d_space,
    const SmartPtr<const SymMatrixSpace> h_space,
    const Matrix& Px_L, const Vector& x_L,
    const Matrix& Px_U, const Vector& x_U,
    Number& df,
    SmartPtr<Vector>& dx,
    SmartPtr<Vector>& dc,
    SmartPtr<Vector>& dd)
  {
    DBG_ASSERT(IsValid(nlp_));
    nlp_->GetScalingParameters(x_space, c_space, d_space, df, dx, dc, dd);
  }

} // namespace Ipopt
