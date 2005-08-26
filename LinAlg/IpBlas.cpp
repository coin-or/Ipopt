// Copyright (C) 2004, 2005 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpBlas.hpp"

// Prototypes for the BLAS routines
extern "C"
{
  /** BLAS Fortran function DDOT */
  double F77_FUNC(ddot,DDOT)(ipfint *n, const double *x, ipfint *incX,
                             const double *y, ipfint *incY);
  /** BLAS Fortran function DNRM2 */
  double F77_FUNC(dnrm2,DNRM2)(ipfint *n, const double *x, ipfint *incX);
  /** BLAS Fortran function DASUM */
  double F77_FUNC(dasum,DASUM)(ipfint *n, const double *x, ipfint *incX);
  /** BLAS Fortran function IDAMAX */
  ipfint F77_FUNC(idamax,IDAMAX)(ipfint *n, const double *x, ipfint *incX);
  /** BLAS Fortran subroutine DCOPY */
  void F77_FUNC(dcopy,DCOPY)(ipfint *n, const double *x, ipfint *incX, double *y,
                             ipfint *incY);
  /** BLAS Fortran subroutine DAXPY */
  void F77_FUNC(daxpy,DAXPY)(ipfint *n, const double *alpha, const double *x,
                             ipfint *incX, double *y, ipfint *incY);
  /** BLAS Fortran subroutine DSCAL */
  void F77_FUNC(dscal,DSCAL)(ipfint *n, const double *alpha, const double *x,
                             ipfint *incX);
}

namespace Ipopt
{
#ifndef HAVE_CBLAS
  /* Interface to FORTRAN routine DDOT. */
  Number IpBlasDdot(Index size, const Number *x, Index incX, const Number *y,
                    Index incY)
  {
    ipfint n=size, INCX=incX, INCY=incY;

    return F77_FUNC(ddot,DDOT)(&n, x, &INCX, y, &INCY);
  }

  /* Interface to FORTRAN routine DNRM2. */
  Number IpBlasDnrm2(Index size, const Number *x, Index incX)
  {
    ipfint n=size, INCX=incX;

    return F77_FUNC(dnrm2,DNRM2)(&n, x, &INCX);
  }

  /* Interface to FORTRAN routine DASUM. */
  Number IpBlasDasum(Index size, const Number *x, Index incX)
  {
    ipfint n=size, INCX=incX;

    return F77_FUNC(dasum,DASUM)(&n, x, &INCX);
  }

  /* Interface to FORTRAN routine DASUM. */
  Index IpBlasIdamax(Index size, const Number *x, Index incX)
  {
    ipfint n=size, INCX=incX;

    return (Index) F77_FUNC(idamax,IDAMAX)(&n, x, &INCX);
  }

  /* Interface to FORTRAN routine DCOPY. */
  void IpBlasDcopy(Index size, const Number *x, Index incX, Number *y, Index incY)
  {
    ipfint N=size, INCX=incX, INCY=incY;

    F77_FUNC(dcopy,DCOPY)(&N, x, &INCX, y, &INCY);
  }

  /* Interface to FORTRAN routine DAXPY. */
  void IpBlasDaxpy(Index size, Number alpha, const Number *x, Index incX, Number *y,
                   Index incY)
  {
    ipfint N=size, INCX=incX, INCY=incY;

    F77_FUNC(daxpy,DAXPY)(&N, &alpha, x, &INCX, y, &INCY);
  }

  /* Interface to FORTRAN routine DSCAL. */
  void IpBlasDscal(Index size, Number alpha, Number *x, Index incX)
  {
    ipfint N=size, INCX=incX;

    F77_FUNC(dscal,DSCAL)(&N, &alpha, x, &INCX);
  }

#else
  /* Interface to CBLAS routine DDOT. */
  Number IpBlasDdot(Index size, const Number *x, Index incX, const Number *y,
                    Index incY)
  {
    Not Implemented Yet!
  }

  /* Interface to CBLAS routine DNRM2. */
  Number IpBlasDnrm2(Index size, const Number *x, Index incX)
  {
    Not Implemented Yet!
  }

  /* Interface to CBLAS routine DASUM. */
  Number IpBlasDasum(Index size, const Number *x, Index incX)
  {
    Not Implemented Yet!
  }

  /* Interface to CBLAS routine DASUM. */
  Index IpBlasIdamax(Index size, const Number *x, Index incX)
  {
    Not Implemented Yet!
  }

  /* Interface to CBLAS routine DCOPY. */
  void IpBlasDcopy(Index size, const Number *x, Index incX, Number *y, Index incY)
  {
    Not Implemented Yet!
  }

  /* Interface to CBLAS routine DAXPY. */
  void IpBlasDaxpy(Index size, Number alpha, const Number *x, Index incX, Number *y,
                   Index incY)
  {
    Not Implemented Yet!
  }

  /* Interface to CBLAS routine DSCAL. */
  void IpBlasDscal(Index size, Number alpha, Number *x, Index incX)
  {
    Not Implemented Yet!
  }
#endif

} // namespace Ipopt
