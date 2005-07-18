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

  void UserScaling::DetermineScalingParametersImpl(
    const SmartPtr<const VectorSpace> x_space,
    const SmartPtr<const VectorSpace> c_space,
    const SmartPtr<const VectorSpace> d_space,
    const SmartPtr<const MatrixSpace> jac_c_space,
    const SmartPtr<const MatrixSpace> jac_d_space,
    const SmartPtr<const SymMatrixSpace> h_space,
    Number& df, Vector& dx,
    Vector& dc, Vector& dd)
  {
    DBG_ASSERT(IsValid(nlp_));
    nlp_->GetScalingParameters(df, dx, dc, dd);
  }

} // namespace Ipopt
