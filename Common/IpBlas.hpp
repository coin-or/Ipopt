// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPBLAS_HPP__
#define __IPBLAS_HPP__

#include "IpUtils.hpp"

namespace Ipopt
{
  // If CBLAS is not available, this is our own interface to the Fortran
  // implementation

  /** Wrapper for BLAS function DDOT.  Compute dot product of vector x
      and vector y */
  Number IpBlasDdot(Index size, const Number *x, Index incX, const Number *y,
                    Index incY);

  /** Wrapper for BLAS function DNRM2.  Compute 2-norm of vector x*/
  Number IpBlasDnrm2(Index size, const Number *x, Index incX);

  /** Wrapper for BLAS function DASUM.  Compute 1-norm of vector x*/
  Number IpBlasDasum(Index size, const Number *x, Index incX);

  /** Wrapper for BLAS function DASUM.  Compute index for largest
      absolute element of vector x */
  Index IpBlasIdamax(Index size, const Number *x, Index incX);

  /** Wrapper for BLAS subroutine DCOPY.  Copying vector x into vector
      y */
  void IpBlasDcopy(Index size, const Number *x, Index incX, Number *y,
                   Index incY);

  /** Wrapper for BLAS subroutine DAXPY.  Adding the alpha multiple of
      vector x to vector y */
  void IpBlasDaxpy(Index size, Number alpha, const Number *x, Index incX,
                   Number *y, Index incY);

  /** Wrapper for BLAS subroutine DSCAL.  Scaling vector x by scalar
      alpha */
  void IpBlasDscal(Index size, Number alpha, Number *x, Index incX);

} // namespace Ipopt

#endif
