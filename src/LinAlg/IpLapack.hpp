// Copyright (C) 2005, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Andreas Waechter              IBM    2005-12-25

#ifndef __IPLAPACK_HPP__
#define __IPLAPACK_HPP__

#include "IpUtils.hpp"
#include "IpException.hpp"

namespace Ipopt
{
DECLARE_STD_EXCEPTION(LAPACK_NOT_INCLUDED);

/** Wrapper for LAPACK subroutine XPOTRS.
 *
 *  Solving a linear system given a Cholesky factorization.
 *  We assume that the Cholesky factor is lower traiangular.
 *  @since 3.14.0
 */
IPOPTLIB_EXPORT void IpLapackPotrs(
   Index         ndim,
   Index         nrhs,
   const Number* a,
   Index         lda,
   Number*       b,
   Index         ldb
);

/** Wrapper for LAPACK subroutine DPOTRS.
 *
 *  Solving a linear system given a Cholesky factorization.
 *  We assume that the Cholesky factor is lower traiangular.
 *
 *  @deprecated Use IpLapackPotrs() instead.
 */
IPOPT_DEPRECATED
inline void IpLapackDpotrs(
   Index         ndim,
   Index         nrhs,
   const Number* a,
   Index         lda,
   Number*       b,
   Index         ldb
)
{
   IpLapackPotrs(ndim, nrhs, a, lda, b, ldb);
}

/** Wrapper for LAPACK subroutine XPOTRF.
 *
 *  Compute Cholesky factorization (lower triangular factor).
 *  info is the return value from the LAPACK routine.
 *  @since 3.14.0
 */
IPOPTLIB_EXPORT void IpLapackPotrf(
   Index   ndim,
   Number* a,
   Index   lda,
   Index&  info
);

/** Wrapper for LAPACK subroutine DPOTRF.
 *
 *  Compute Cholesky factorization (lower triangular factor).
 *  info is the return value from the LAPACK routine.
 *
 *  @deprecated Use IpLapackPotrf() instead.
 */
IPOPT_DEPRECATED
inline void IpLapackDpotrf(
   Index   ndim,
   Number* a,
   Index   lda,
   Index&  info
)
{
   IpLapackPotrf(ndim, a, lda, info);
}

/** Wrapper for LAPACK subroutine XSYEV.
 *
 *  Compute the Eigenvalue decomposition for a given matrix.
 *  If compute_eigenvectors is true, a will contain the eigenvectors
 *  in its columns on return.
 *  @since 3.14.0
 */
IPOPTLIB_EXPORT void IpLapackSyev(
   bool    compute_eigenvectors,
   Index   ndim,
   Number* a,
   Index   lda,
   Number* w,
   Index&  info
);

/** Wrapper for LAPACK subroutine DSYEV.
 *
 *  Compute the Eigenvalue decomposition for a given matrix.
 *  If compute_eigenvectors is true, a will contain the eigenvectors
 *  in its columns on return.
 *
 *  @deprecated Use IpLapackSyev() instead
 */
IPOPT_DEPRECATED
inline void IpLapackDsyev(
   bool    compute_eigenvectors,
   Index   ndim,
   Number* a,
   Index   lda,
   Number* w,
   Index&  info
)
{
   IpLapackSyev(compute_eigenvectors, ndim, a, lda, w, info);
}

/** Wrapper for LAPACK subroutine XGETRF.
 *
 *  Compute LU factorization.
 *  info is the return value from the LAPACK routine.
 *  @since 3.14.0
 */
IPOPTLIB_EXPORT void IpLapackGetrf(
   Index   ndim,
   Number* a,
   Index*  ipiv,
   Index   lda,
   Index&  info
);

/** Wrapper for LAPACK subroutine DGETRF.
 *
 *  Compute LU factorization.
 *  info is the return value from the LAPACK routine.
 *
 *  @deprecated Use IpLapackGetrf() instead.
 */
IPOPT_DEPRECATED
inline void IpLapackDgetrf(
   Index   ndim,
   Number* a,
   Index*  ipiv,
   Index   lda,
   Index&  info
)
{
   IpLapackGetrf(ndim, a, ipiv, lda, info);
}

/** Wrapper for LAPACK subroutine XGETRS.
 *
 * Solving a linear system given a LU factorization.
 * @since 3.14.0
 */
IPOPTLIB_EXPORT void IpLapackGetrs(
   Index         ndim,
   Index         nrhs,
   const Number* a,
   Index         lda,
   Index*        ipiv,
   Number*       b,
   Index         ldb
);

/** Wrapper for LAPACK subroutine DGETRS.
 *
 * Solving a linear system given a LU factorization.
 *
 * @deprecated Use IpLapackGetrs() instead.
 */
IPOPT_DEPRECATED
inline void IpLapackDgetrs(
   Index         ndim,
   Index         nrhs,
   const Number* a,
   Index         lda,
   Index*        ipiv,
   Number*       b,
   Index         ldb
)
{
   IpLapackGetrs(ndim, nrhs, a, lda, ipiv, b, ldb);
}

/** Wrapper for LAPACK subroutine XPPSV.
 *
 *  Solves a symmetric positive
 *  definite linear system in packed storage format (upper triangular).
 *  info is the return value from the LAPACK routine.
 *  @since 3.14.0
 */
IPOPTLIB_EXPORT void IpLapackPpsv(
   Index         ndim,
   Index         nrhs,
   const Number* a,
   Number*       b,
   Index         ldb,
   Index&        info
);

/** Wrapper for LAPACK subroutine DPPSV.
 *
 *  Solves a symmetric positive
 *  definite linear system in packed storage format (upper triangular).
 *  info is the return value from the LAPACK routine.
 *
 *  @deprecated Use IpLapackPpsv() instead.
 */
IPOPT_DEPRECATED
inline void IpLapackDppsv(
   Index         ndim,
   Index         nrhs,
   const Number* a,
   Number*       b,
   Index         ldb,
   Index&        info
)
{
   IpLapackPpsv(ndim, nrhs, a, b, ldb, info);
}

} // namespace Ipopt

#endif
