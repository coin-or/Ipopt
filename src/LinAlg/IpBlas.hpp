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
/** Wrapper for BLAS function XDOT.
 *
 * Compute dot product of vector x and vector y.
 */
IPOPTLIB_EXPORT Number IpBlasDot(
   Index         size,
   const Number* x,
   Index         incX,
   const Number* y,
   Index         incY
);

/** Wrapper for BLAS function XNRM2.
 *
 * Compute 2-norm of vector x.
 */
IPOPTLIB_EXPORT Number IpBlasNrm2(
   Index         size,
   const Number* x,
   Index         incX
);

/** Wrapper for BLAS function XASUM.
 *
 * Compute 1-norm of vector x.
 */
IPOPTLIB_EXPORT Number IpBlasAsum(
   Index         size,
   const Number* x,
   Index         incX
);

/** Wrapper for BLAS function IXAMAX.
 *
 * Compute index for largest absolute element of vector x.
 */
IPOPTLIB_EXPORT int IpBlasIamax(
   Index         size,
   const Number* x,
   Index         incX
);

/** Wrapper for BLAS subroutine XCOPY.
 *
 * Copying vector x into vector y.
 */
IPOPTLIB_EXPORT void IpBlasCopy(
   Index         size,
   const Number* x,
   Index         incX,
   Number*       y,
   Index         incY
);

/** Wrapper for BLAS subroutine XAXPY.
 *
 * Adding the alpha multiple of vector x to vector y.
 */
IPOPTLIB_EXPORT void IpBlasAxpy(
   Index         size,
   Number        alpha,
   const Number* x,
   Index         incX,
   Number*       y,
   Index         incY
);

/** Wrapper for BLAS subroutine XSCAL.
 *
 * Scaling vector x by scalar alpha.
 */
IPOPTLIB_EXPORT void IpBlasScal(
   Index   size,
   Number  alpha,
   Number* x,
   Index   incX
);

/** Wrapper for BLAS subroutine XGEMV.
 *
 * Multiplying a matrix with a vector.
 */
IPOPTLIB_EXPORT void IpBlasGemv(
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

/** Wrapper for BLAS subroutine XSYMV.
 *
 * Multiplying a symmetric matrix with a vector.
 */
IPOPTLIB_EXPORT void IpBlasSymv(
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

/** Wrapper for BLAS subroutine XGEMM.
 *
 * Multiplying two matrices.
 */
IPOPTLIB_EXPORT void IpBlasGemm(
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

/** Wrapper for BLAS subroutine XSYRK.
 *
 * Adding a high-rank update to a matrix.
 */
IPOPTLIB_EXPORT void IpBlasSyrk(
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

/** Wrapper for BLAS subroutine XTRSM.
 *
 * Backsolve for a lower triangular matrix.
 */
IPOPTLIB_EXPORT void IpBlasTrsm(
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
