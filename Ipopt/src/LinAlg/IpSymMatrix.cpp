// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpSymMatrix.hpp"

namespace Ipopt
{
  void SymMatrix::TransMultVectorImpl(Number alpha, const Vector& x, Number beta,
                                      Vector& y) const
  {
    // Since this matrix is symetric, this is the same operation as
    // MultVector
    MultVector(alpha, x, beta, y);
  }
} // namespace Ipopt
