// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpoptConfig.h"
#include "IpBlas.hpp"

/* we currently have no separate check for Blas, but assume that Blas comes with Lapack
 * thus, we use the nameing convention of Lapack for Blas, too
 */
#ifndef IPOPT_BLAS_FUNC
#define IPOPT_BLAS_FUNC(name,NAME) IPOPT_LAPACK_FUNC(name,NAME)
#endif

#include <cstring>

// Prototypes for the BLAS routines
extern "C"
{
#ifdef IPOPT_SINGLE
   /** BLAS Fortran function SDOT */
   float IPOPT_BLAS_FUNC(sdot, SDOT)(
      ipfint*      n,
      const float* x,
      ipfint*      incX,
      const float* y,
      ipfint*      incY
   );

   /** BLAS Fortran function SNRM2 */
   float IPOPT_BLAS_FUNC(snrm2, SNRM2)(
      ipfint*      n,
      const float* x,
      ipfint*      incX
   );

   /** BLAS Fortran function SASUM */
   float IPOPT_BLAS_FUNC(sasum, SASUM)(
      ipfint*      n,
      const float* x,
      ipfint*      incX
   );

   /** BLAS Fortran function ISAMAX */
   ipfint IPOPT_BLAS_FUNC(isamax, ISAMAX)(
      ipfint*      n,
      const float* x,
      ipfint*      incX
   );

   /** BLAS Fortran subroutine SCOPY */
   void IPOPT_BLAS_FUNC(scopy, SCOPY)(
      ipfint*      n,
      const float* x,
      ipfint*      incX,
      float*       y,
      ipfint*      incY
   );

   /** BLAS Fortran subroutine SAXPY */
   void IPOPT_BLAS_FUNC(saxpy, SAXPY)(
      ipfint*      n,
      const float* alpha,
      const float* x,
      ipfint*      incX,
      float*       y,
      ipfint*      incY
   );

   /** BLAS Fortran subroutine SSCAL */
   void IPOPT_BLAS_FUNC(sscal, SSCAL)(
      ipfint*      n,
      const float* alpha,
      const float* x,
      ipfint*      incX
   );

   /** BLAS Fortran subroutine SGEMV */
   void IPOPT_BLAS_FUNC(sgemv, SGEMV)(
      char*        trans,
      ipfint*      m,
      ipfint*      n,
      const float* alpha,
      const float* a,
      ipfint*      lda,
      const float* x,
      ipfint*      incX,
      const float* beta,
      float*       y,
      ipfint*      incY,
      int          trans_len
   );

   /** BLAS Fortran subroutine SSYMV */
   void IPOPT_BLAS_FUNC(ssymv, SSYMV)(
      char*        uplo,
      ipfint*      n,
      const float* alpha,
      const float* a,
      ipfint*      lda,
      const float* x,
      ipfint*      incX,
      const float* beta,
      float*       y,
      ipfint*      incY,
      int          uplo_len
   );

   /** BLAS Fortran subroutine SGEMM */
   void IPOPT_BLAS_FUNC(sgemm, SGEMM)(
      char*        transa,
      char*        transb,
      ipfint*      m,
      ipfint*      n,
      ipfint*      k,
      const float* alpha,
      const float* a,
      ipfint*      lda,
      const float* b,
      ipfint*      ldb,
      const float* beta,
      float*       c,
      ipfint*      ldc,
      int          transa_len,
      int          transb_len
   );

   /** BLAS Fortran subroutine SSYRK */
   void IPOPT_BLAS_FUNC(ssyrk, SSYRK)(
      char*        uplo,
      char*        trans,
      ipfint*      n,
      ipfint*      k,
      const float* alpha,
      const float* a,
      ipfint*      lda,
      const float* beta,
      float*       c,
      ipfint*      ldc,
      int          uplo_len,
      int          trans_len
   );

   /** BLAS Fortran subroutine STRSM */
   void IPOPT_BLAS_FUNC(strsm, STRSM)(
      char*        side,
      char*        uplo,
      char*        transa,
      char*        diag,
      ipfint*      m,
      ipfint*      n,
      const float* alpha,
      const float* a,
      ipfint*      lda,
      const float* b,
      ipfint*      ldb,
      int          side_len,
      int          uplo_len,
      int          transa_len,
      int          diag_len
   );
#else
   /** BLAS Fortran function DDOT */
   double IPOPT_BLAS_FUNC(ddot, DDOT)(
      ipfint*       n,
      const double* x,
      ipfint*       incX,
      const double* y,
      ipfint*       incY
   );

   /** BLAS Fortran function DNRM2 */
   double IPOPT_BLAS_FUNC(dnrm2, DNRM2)(
      ipfint*       n,
      const double* x,
      ipfint*       incX
   );

   /** BLAS Fortran function DASUM */
   double IPOPT_BLAS_FUNC(dasum, DASUM)(
      ipfint*       n,
      const double* x,
      ipfint*       incX
   );

   /** BLAS Fortran function IDAMAX */
   ipfint IPOPT_BLAS_FUNC(idamax, IDAMAX)(
      ipfint*       n,
      const double* x,
      ipfint*       incX
   );

   /** BLAS Fortran subroutine DCOPY */
   void IPOPT_BLAS_FUNC(dcopy, DCOPY)(
      ipfint*       n,
      const double* x,
      ipfint*       incX,
      double*       y,
      ipfint*       incY
   );

   /** BLAS Fortran subroutine DAXPY */
   void IPOPT_BLAS_FUNC(daxpy, DAXPY)(
      ipfint*       n,
      const double* alpha,
      const double* x,
      ipfint*       incX,
      double*       y,
      ipfint*       incY
   );

   /** BLAS Fortran subroutine DSCAL */
   void IPOPT_BLAS_FUNC(dscal, DSCAL)(
      ipfint*       n,
      const double* alpha,
      const double* x,
      ipfint*       incX
   );

   /** BLAS Fortran subroutine DGEMV */
   void IPOPT_BLAS_FUNC(dgemv, DGEMV)(
      char*         trans,
      ipfint*       m,
      ipfint*       n,
      const double* alpha,
      const double* a,
      ipfint*       lda,
      const double* x,
      ipfint*       incX,
      const double* beta,
      double*       y,
      ipfint*       incY,
      int           trans_len
   );

   /** BLAS Fortran subroutine DSYMV */
   void IPOPT_BLAS_FUNC(dsymv, DSYMV)(
      char*         uplo,
      ipfint*       n,
      const double* alpha,
      const double* a,
      ipfint*       lda,
      const double* x,
      ipfint*       incX,
      const double* beta,
      double*       y,
      ipfint*       incY,
      int           uplo_len
   );

   /** BLAS Fortran subroutine DGEMM */
   void IPOPT_BLAS_FUNC(dgemm, DGEMM)(
      char*         transa,
      char*         transb,
      ipfint*       m,
      ipfint*       n,
      ipfint*       k,
      const double* alpha,
      const double* a,
      ipfint*       lda,
      const double* b,
      ipfint*       ldb,
      const double* beta,
      double*       c,
      ipfint*       ldc,
      int           transa_len,
      int           transb_len
   );

   /** BLAS Fortran subroutine DSYRK */
   void IPOPT_BLAS_FUNC(dsyrk, DSYRK)(
      char*         uplo,
      char*         trans,
      ipfint*       n,
      ipfint*       k,
      const double* alpha,
      const double* a,
      ipfint*       lda,
      const double* beta,
      double*       c,
      ipfint*       ldc,
      int           uplo_len,
      int           trans_len
   );

   /** BLAS Fortran subroutine DTRSM */
   void IPOPT_BLAS_FUNC(dtrsm, DTRSM)(
      char*         side,
      char*         uplo,
      char*         transa,
      char*         diag,
      ipfint*       m,
      ipfint*       n,
      const double* alpha,
      const double* a,
      ipfint*       lda,
      const double* b,
      ipfint*       ldb,
      int           side_len,
      int           uplo_len,
      int           transa_len,
      int           diag_len
   );

#endif
}
namespace Ipopt
{
Number IpBlasDot(
   Index         size,
   const Number* x,
   Index         incX,
   const Number* y,
   Index         incY
)
{
   if( incX > 0 && incY > 0 )
   {
      ipfint n = size, INCX = incX, INCY = incY;
#ifdef IPOPT_SINGLE
      return IPOPT_BLAS_FUNC(sdot, SDOT)(&n, x, &INCX, y, &INCY);
#else
      return IPOPT_BLAS_FUNC(ddot, DDOT)(&n, x, &INCX, y, &INCY);
#endif
   }
   else
   {
      Number s = 0.0;

      for( ; size; --size, x += incX, y += incY )
      {
         s += *x * *y;
      }

      return s;
   }
}

Number IpBlasNrm2(
   Index         size,
   const Number* x,
   Index         incX
)
{
   ipfint n = size, INCX = incX;
#ifdef IPOPT_SINGLE
   return IPOPT_BLAS_FUNC(snrm2, SNRM2)(&n, x, &INCX);
#else
   return IPOPT_BLAS_FUNC(dnrm2, DNRM2)(&n, x, &INCX);
#endif
}

Number IpBlasAsum(
   Index         size,
   const Number* x,
   Index         incX
)
{
   ipfint n = size, INCX = incX;
#ifdef IPOPT_SINGLE
   return IPOPT_BLAS_FUNC(sasum, SASUM)(&n, x, &INCX);
#else
   return IPOPT_BLAS_FUNC(dasum, DASUM)(&n, x, &INCX);
#endif
}

/** interface to FORTRAN routine IXAMAX */
Index IpBlasIamax(
   Index         size,
   const Number* x,
   Index         incX
)
{
   ipfint n = size, INCX = incX;
#ifdef IPOPT_SINGLE
   return (Index) IPOPT_BLAS_FUNC(isamax, ISAMAX)(&n, x, &INCX);
#else
   return (Index) IPOPT_BLAS_FUNC(idamax, IDAMAX)(&n, x, &INCX);
#endif
}

/** interface to FORTRAN routine XCOPY */
void IpBlasCopy(
   Index         size,
   const Number* x,
   Index         incX,
   Number*       y,
   Index         incY
)
{
   if( incX > 0 )
   {
      ipfint N = size, INCX = incX, INCY = incY;
#ifdef IPOPT_SINGLE
      IPOPT_BLAS_FUNC(scopy, SCOPY)(&N, x, &INCX, y, &INCY);
#else
      IPOPT_BLAS_FUNC(dcopy, DCOPY)(&N, x, &INCX, y, &INCY);
#endif
   }
   else if( incY == 1 )
   {
      for( ; size; --size, ++y )
      {
         *y = *x;
      }
   }
   else
   {
      for( ; size; --size, y += incY )
      {
         *y = *x;
      }
   }
}

void IpBlasAxpy(
   Index         size,
   Number        alpha,
   const Number* x,
   Index         incX,
   Number*       y,
   Index         incY
)
{
   if( incX > 0 )
   {
      ipfint N = size, INCX = incX, INCY = incY;
#ifdef IPOPT_SINGLE
      IPOPT_BLAS_FUNC(saxpy, SAXPY)(&N, &alpha, x, &INCX, y, &INCY);
#else
      IPOPT_BLAS_FUNC(daxpy, DAXPY)(&N, &alpha, x, &INCX, y, &INCY);
#endif
   }
   else if( incY == 1 )
   {
      for( ; size; --size, ++y )
      {
         *y += alpha * *x;
      }
   }
   else
   {
      for( ; size; --size, y += incY )
      {
         *y += alpha * *x;
      }
   }
}

void IpBlasScal(
   Index   size,
   Number  alpha,
   Number* x,
   Index   incX
)
{
   ipfint N = size, INCX = incX;
#ifdef IPOPT_SINGLE
   IPOPT_BLAS_FUNC(sscal, SSCAL)(&N, &alpha, x, &INCX);
#else
   IPOPT_BLAS_FUNC(dscal, DSCAL)(&N, &alpha, x, &INCX);
#endif
}

void IpBlasGemv(
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
)
{
   ipfint M = nCols, N = nRows, LDA = ldA, INCX = incX, INCY = incY;

   char TRANS;
   if( trans )
   {
      TRANS = 'T';
   }
   else
   {
      TRANS = 'N';
   }
#ifdef IPOPT_SINGLE
   IPOPT_BLAS_FUNC(sgemv, SGEMV)(&TRANS, &M, &N, &alpha, A, &LDA, x, &INCX, &beta, y, &INCY, 1);
#else
   IPOPT_BLAS_FUNC(dgemv, DGEMV)(&TRANS, &M, &N, &alpha, A, &LDA, x, &INCX, &beta, y, &INCY, 1);
#endif
}

void IpBlasSymv(
   Index         n,
   Number        alpha,
   const Number* A,
   Index         ldA,
   const Number* x,
   Index         incX,
   Number        beta,
   Number*       y,
   Index         incY
)
{
   ipfint N = n, LDA = ldA, INCX = incX, INCY = incY;

   char UPLO = 'L';
#ifdef IPOPT_SINGLE
   IPOPT_BLAS_FUNC(ssymv, SSYMV)(&UPLO, &N, &alpha, A, &LDA, x, &INCX, &beta, y, &INCY, 1);
#else
   IPOPT_BLAS_FUNC(dsymv, DSYMV)(&UPLO, &N, &alpha, A, &LDA, x, &INCX, &beta, y, &INCY, 1);
#endif
}

void IpBlasGemm(
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
)
{
   ipfint M = m, N = n, K = k, LDA = ldA, LDB = ldB, LDC = ldC;

   char TRANSA;
   if( transa )
   {
      TRANSA = 'T';
   }
   else
   {
      TRANSA = 'N';
   }
   char TRANSB;
   if( transb )
   {
      TRANSB = 'T';
   }
   else
   {
      TRANSB = 'N';
   }
#ifdef IPOPT_SINGLE
   IPOPT_BLAS_FUNC(sgemm, SGEMM)(&TRANSA, &TRANSB, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC, 1, 1);
#else
   IPOPT_BLAS_FUNC(dgemm, DGEMM)(&TRANSA, &TRANSB, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC, 1, 1);
#endif
}

void IpBlasSyrk(
   bool          trans,
   Index         ndim,
   Index         nrank,
   Number        alpha,
   const Number* A,
   Index         ldA,
   Number        beta,
   Number*       C,
   Index         ldC
)
{
   ipfint N = ndim, K = nrank, LDA = ldA, LDC = ldC;

   char UPLO = 'L';
   char TRANS;
   if( trans )
   {
      TRANS = 'T';
   }
   else
   {
      TRANS = 'N';
   }
#ifdef IPOPT_SINGLE
   IPOPT_BLAS_FUNC(ssyrk, SSYRK)(&UPLO, &TRANS, &N, &K, &alpha, A, &LDA, &beta, C, &LDC, 1, 1);
#else
   IPOPT_BLAS_FUNC(dsyrk, DSYRK)(&UPLO, &TRANS, &N, &K, &alpha, A, &LDA, &beta, C, &LDC, 1, 1);
#endif
}

void IpBlasTrsm(
   bool          trans,
   Index         ndim,
   Index         nrhs,
   Number        alpha,
   const Number* A,
   Index         ldA,
   Number*       B,
   Index         ldB
)
{
   ipfint M = ndim, N = nrhs, LDA = ldA, LDB = ldB;

   char SIDE = 'L';
   char UPLO = 'L';
   char TRANSA;
   if( trans )
   {
      TRANSA = 'T';
   }
   else
   {
      TRANSA = 'N';
   }
   char DIAG = 'N';
#ifdef IPOPT_SINGLE
   IPOPT_BLAS_FUNC(strsm, STRSM)(&SIDE, &UPLO, &TRANSA, &DIAG, &M, &N, &alpha, A, &LDA, B, &LDB, 1, 1, 1, 1);
#else
   IPOPT_BLAS_FUNC(dtrsm, DTRSM)(&SIDE, &UPLO, &TRANSA, &DIAG, &M, &N, &alpha, A, &LDA, B, &LDB, 1, 1, 1, 1);
#endif
}

} // namespace Ipopt
