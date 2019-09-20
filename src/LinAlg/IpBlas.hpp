// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPBLAS_HPP__
#define __IPBLAS_HPP__

#include "IpUtils.hpp"

namespace Ipopt
{
/** Wrapper for BLAS function DDOT.
 *
 * Compute dot product of vector x and vector y.
 */
IPOPTLIB_EXPORT Number IpBlasDdot(
   Index         size,
   const Number* x,
   Index         incX,
   const Number* y,
   Index         incY
);

/** Wrapper for BLAS function DNRM2.
 *
 * Compute 2-norm of vector x.
 */
IPOPTLIB_EXPORT Number IpBlasDnrm2(
   Index         size,
   const Number* x,
   Index         incX
);

/** Wrapper for BLAS function DASUM.
 *
 * Compute 1-norm of vector x.
 */
IPOPTLIB_EXPORT Number IpBlasDasum(
   Index         size,
   const Number* x,
   Index         incX
);

/** Wrapper for BLAS function IDAMAX.
 *
 * Compute index for largest absolute element of vector x.
 */
IPOPTLIB_EXPORT int IpBlasIdamax(
   Index         size,
   const Number* x,
   Index         incX
);

/** Wrapper for BLAS subroutine DCOPY.
 *
 * Copying vector x into vector y.
 */
IPOPTLIB_EXPORT void IpBlasDcopy(
   Index         size,
   const Number* x,
   Index         incX,
   Number*       y,
   Index         incY
);

/** Wrapper for BLAS subroutine DAXPY.
 *
 * Adding the alpha multiple of vector x to vector y.
 */
IPOPTLIB_EXPORT void IpBlasDaxpy(
   Index         size,
   Number        alpha,
   const Number* x,
   Index         incX,
   Number*       y,
   Index         incY
);

/** Wrapper for BLAS subroutine DSCAL.
 *
 * Scaling vector x by scalar alpha.
 */
IPOPTLIB_EXPORT void IpBlasDscal(
   Index   size,
   Number  alpha,
   Number* x,
   Index   incX
);

/** Wrapper for BLAS subroutine DGEMV.
 *
 * Multiplying a matrix with a vector.
 */
IPOPTLIB_EXPORT void IpBlasDgemv(
   bool          trans,
   Index         nRows,
   Index         nCols,
   Number        alpha,
   const Number* A,
   Index         ldA,
   const Number* x,
   Index         incX,
   Number        beta,
   Number*       y,
   Index         incY
);

/** Wrapper for BLAS subroutine DSYMV.
 *
 * Multiplying a symmetric matrix with a vector.
 */
IPOPTLIB_EXPORT void IpBlasDsymv(
   Index         n,
   Number        alpha,
   const Number* A,
   Index         ldA,
   const Number* x,
   Index         incX,
   Number        beta,
   Number*       y,
   Index         incY
);

/** Wrapper for BLAS subroutine DGEMM.
 *
 * Multiplying two matrices.
 */
IPOPTLIB_EXPORT void IpBlasDgemm(
   bool          transa,
   bool          transb,
   Index         m,
   Index         n,
   Index         k,
   Number        alpha,
   const Number* A,
   Index         ldA,
   const Number* B,
   Index         ldB,
   Number        beta,
   Number*       C,
   Index         ldC
);

/** Wrapper for BLAS subroutine DSYRK.
 *
 * Adding a high-rank update to a matrix.
 */
IPOPTLIB_EXPORT void IpBlasDsyrk(
   bool          trans,
   Index         ndim,
   Index         nrank,
   Number        alpha,
   const Number* A,
   Index         ldA,
   Number        beta,
   Number*       C,
   Index         ldC
);

/** Wrapper for BLAS subroutine DTRSM.
 *
 * Backsolve for a lower triangular matrix.
 */
IPOPTLIB_EXPORT void IpBlasDtrsm(
   bool          trans,
   Index         ndim,
   Index         nrhs,
   Number        alpha,
   const Number* A,
   Index         ldA,
   Number*       B,
   Index         ldB
);

} // namespace Ipopt

#endif
