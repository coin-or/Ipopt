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
                                 Number b, const Vector& v2, Number c)
  {
    if (c==0.) {
      if (a==1.) {
        Copy(v1);
        if (b!=0.) {
          Axpy(b, v2);
        }
      }
      else if (a==0.) {
        if (b==0.) {
          Set(0.);
        }
        else {
          Copy(v2);
          if (b!=1.) {
            Scal(b);
          }
        }
      }
      else {
        if (b==1.) {
          Copy(v2);
          Axpy(a, v1);
        }
        else if (b==0.) {
          Copy(v1);
          Scal(a);
        }
        else {
          Copy(v1);
          Scal(a);
          Axpy(b, v2);
        }
      }
    }
    else { /* c==0. */
      if (c!=1.) {
        Scal(c);
      }
      if (a!=0.) {
        Axpy(a, v1);
      }
      if (b!=0.) {
        Axpy(b, v2);
      }
    }
  }

} // namespace Ipopt
