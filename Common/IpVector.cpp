// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpVector.hpp"

namespace Ipopt
{

  Vector::Vector(const VectorSpace* owner_space)
      :
      owner_space_(owner_space),
      TaggedObject(),
      dot_cache_(10),
      nrm2_cache_(1),
      asum_cache_(1),
      amax_cache_(1),
      sum_cache_(1),
      sumlogs_cache_(1)
  {
    DBG_ASSERT(IsValid(owner_space_));
  }

  Vector::~Vector()
  {}

  VectorSpace::VectorSpace(Index dim)
      :
      dim_(dim)
  {}

  VectorSpace::~VectorSpace()
  {}

  /* Prototype implementation for specialized functions */
  void Vector::AddTwoVectorsImpl(Number a, const Vector& v1,
				 Number b, const Vector& v2)
  {
    Copy(v1);
    if (a!=1.) {
      Scal(a);
    }
    if (b!=0.) {
      Axpy(b, v2);
    }
  }

  /*
  Vector* VectorSpace::MakeNewEVector_Vector(Number factor)
  {
    return MakeNewEVector(factor);
  }

  EVector* VectorSpace::MakeNewEVector(Number factor)
  {
    EVector* retVector = new EVector(this);
    retVector->Set(factor);
    return retVector;
  }
  */

} // namespace Ipopt
