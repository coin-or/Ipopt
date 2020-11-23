// Copyright (C) 2005, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Andreas Waechter              IBM    2005-12-25

#include "IpoptConfig.h"
#include "IpLapack.hpp"

#ifdef FUNNY_LAPACK_FINT
# define ipfint long
# define ipfintarray int
#else
# define ipfintarray ipfint
#endif

#ifdef IPOPT_HAS_LAPACK
// Prototypes for the LAPACK routines
#ifdef IPOPT_SINGLE
extern "C"
{
   /** LAPACK Fortran subroutine SPOTRS. */
   void IPOPT_LAPACK_FUNC(spotrs, SPOTRS)(
      char*        uplo,
      ipfint*      n,
      ipfint*      nrhs,
      const float* A,
      ipfint*      ldA,
      float*       B,
      ipfint*      ldB,
      ipfint*      info,
      int          uplo_len
   );

   /** LAPACK Fortran subroutine SPOTRF. */
   void IPOPT_LAPACK_FUNC(spotrf, SPOTRF)(
      char*   uplo,
      ipfint* n,
      float*  A,
      ipfint* ldA,
      ipfint* info,
      int     uplo_len
   );

   /** LAPACK Fortran subroutine SSYEV */
   void IPOPT_LAPACK_FUNC(ssyev, SSYEV)(
      char*   jobz,
      char*   uplo,
      ipfint* n,
      float*  A,
      ipfint* ldA,
      float*  W,
      float*  WORK,
      ipfint* LWORK,
      ipfint* info,
      int     jobz_len,
      int     uplo_len
   );

   /** LAPACK Fortran subroutine SGETRF. */
   void IPOPT_LAPACK_FUNC(sgetrf, SGETRF)(
      ipfint*      m,
      ipfint*      n,
      float*       A,
      ipfint*      ldA,
      ipfintarray* IPIV,
      ipfint*      info
   );

   /** LAPACK Fortran subroutine SGETRS. */
   void IPOPT_LAPACK_FUNC(sgetrs, SGETRS)(
      char*        trans,
      ipfint*      n,
      ipfint*      nrhs,
      const float* A,
      ipfint*      ldA,
      ipfintarray* IPIV,
      float*       B,
      ipfint*      ldB,
      ipfint*      info,
      int          trans_len
   );

   /** LAPACK Fortran subroutine SPPSV. */
   void IPOPT_LAPACK_FUNC(sppsv, SPPSV)(
      char*        uplo,
      ipfint*      n,
      ipfint*      nrhs,
      const float* A,
      float*       B,
      ipfint*      ldB,
      ipfint*      info
   );
}

#else
extern "C"
{
   /** LAPACK Fortran subroutine DPOTRS. */
   void IPOPT_LAPACK_FUNC(dpotrs, DPOTRS)(
      char*         uplo,
      ipfint*       n,
      ipfint*       nrhs,
      const double* A,
      ipfint*       ldA,
      double*       B,
      ipfint*       ldB,
      ipfint*       info,
      int           uplo_len
   );

   /** LAPACK Fortran subroutine DPOTRF. */
   void IPOPT_LAPACK_FUNC(dpotrf, DPOTRF)(
      char*   uplo,
      ipfint* n,
      double* A,
      ipfint* ldA,
      ipfint* info,
      int     uplo_len
   );

   /** LAPACK Fortran subroutine DSYEV */
   void IPOPT_LAPACK_FUNC(dsyev, DSYEV)(
      char*   jobz,
      char*   uplo,
      ipfint* n,
      double* A,
      ipfint* ldA,
      double* W,
      double* WORK,
      ipfint* LWORK,
      ipfint* info,
      int     jobz_len,
      int     uplo_len
   );

   /** LAPACK Fortran subroutine DGETRF. */
   void IPOPT_LAPACK_FUNC(dgetrf, DGETRF)(
      ipfint*      m,
      ipfint*      n,
      double*      A,
      ipfint*      ldA,
      ipfintarray* IPIV,
      ipfint*      info
   );

   /** LAPACK Fortran subroutine DGETRS. */
   void IPOPT_LAPACK_FUNC(dgetrs, DGETRS)(
      char*         trans,
      ipfint*       n,
      ipfint*       nrhs,
      const double* A,
      ipfint*       ldA,
      ipfintarray*  IPIV,
      double*       B,
      ipfint*       ldB,
      ipfint*       info,
      int           trans_len
   );

   /** LAPACK Fortran subroutine DPPSV. */
   void IPOPT_LAPACK_FUNC(dppsv, DPPSV)(
      char*         uplo,
      ipfint*       n,
      ipfint*       nrhs,
      const double* A,
      double*       B,
      ipfint*       ldB,
      ipfint*       info
   );
}
#endif  /* ifdef IPOPT_SINGLE */

namespace Ipopt
{
void IpLapackPotrs(
   Index         ndim,
   Index         nrhs,
   const Number* a,
   Index         lda,
   Number*       b,
   Index         ldb
)
{
#ifdef IPOPT_HAS_LAPACK
   ipfint N = ndim, NRHS = nrhs, LDA = lda, LDB = ldb, INFO;
   char uplo = 'L';

#ifdef IPOPT_SINGLE
   IPOPT_LAPACK_FUNC(spotrs, SPOTRS)(&uplo, &N, &NRHS, a, &LDA, b, &LDB, &INFO, 1);
#else
   IPOPT_LAPACK_FUNC(dpotrs, DPOTRS)(&uplo, &N, &NRHS, a, &LDA, b, &LDB, &INFO, 1);
#endif
   DBG_ASSERT(INFO == 0);
#else

   std::string msg =
      "Ipopt has been compiled without LAPACK routine DPOTRS, but options are chosen that require this dependency.  Abort.";
   THROW_EXCEPTION(LAPACK_NOT_INCLUDED, msg);
#endif
}

void IpLapackPotrf(
   Index   ndim,
   Number* a,
   Index   lda,
   Index&  info
)
{
#ifdef IPOPT_HAS_LAPACK
   ipfint N = ndim, LDA = lda, INFO;

   char UPLO = 'L';

#ifdef IPOPT_SINGLE
   IPOPT_LAPACK_FUNC(spotrf, SPOTRF)(&UPLO, &N, a, &LDA, &INFO, 1);
#else
   IPOPT_LAPACK_FUNC(dpotrf, DPOTRF)(&UPLO, &N, a, &LDA, &INFO, 1);
#endif

   info = INFO;
#else

   std::string msg =
      "Ipopt has been compiled without LAPACK routine DPOTRF, but options are chosen that require this dependency.  Abort.";
   THROW_EXCEPTION(LAPACK_NOT_INCLUDED, msg);
#endif

}

void IpLapackSyev(
   bool    compute_eigenvectors,
   Index   ndim,
   Number* a,
   Index   lda,
   Number* w,
   Index&  info
)
{
#ifdef IPOPT_HAS_LAPACK
   ipfint N = ndim, LDA = lda, INFO;

   char JOBZ;
   if (compute_eigenvectors)
   {
      JOBZ = 'V';
   }
   else
   {
      JOBZ = 'N';
   }
   char UPLO = 'L';

   // First we find out how large LWORK should be
   ipfint LWORK = -1;
   Number WORK_PROBE;
#ifdef IPOPT_SINGLE
   IPOPT_LAPACK_FUNC(ssyev, SSYEV)(&JOBZ, &UPLO, &N, a, &LDA, w,
                                  &WORK_PROBE, &LWORK, &INFO, 1, 1);
#else
   IPOPT_LAPACK_FUNC(dsyev, DSYEV)(&JOBZ, &UPLO, &N, a, &LDA, w,
                                  &WORK_PROBE, &LWORK, &INFO, 1, 1);
#endif
   DBG_ASSERT(INFO == 0);

   LWORK = (ipfint) WORK_PROBE;
   DBG_ASSERT(LWORK > 0);

   Number* WORK = new Number[LWORK];
   for (Index i = 0; i < LWORK; i++)
   {
      WORK[i] = i;
   }
#ifdef IPOPT_SINGLE
   IPOPT_LAPACK_FUNC(ssyev, SSYEV)(&JOBZ, &UPLO, &N, a, &LDA, w,
                                  WORK, &LWORK, &INFO, 1, 1);
#else
   IPOPT_LAPACK_FUNC(dsyev, DSYEV)(&JOBZ, &UPLO, &N, a, &LDA, w,
                                  WORK, &LWORK, &INFO, 1, 1);
#endif

   DBG_ASSERT(INFO >= 0);
   info = INFO;

   delete [] WORK;
#else

   std::string msg =
      "Ipopt has been compiled without LAPACK routine DSYEV, but options are chosen that require this dependency.  Abort.";
   THROW_EXCEPTION(LAPACK_NOT_INCLUDED, msg);
#endif

}

void IpLapackGetrf(
   Index   ndim,
   Number* a,
   Index*  ipiv,
   Index   lda,
   Index&  info
)
{
#ifdef IPOPT_HAS_LAPACK
   ipfint M = ndim, N = ndim, LDA = lda, INFO;

#ifdef IPOPT_SINGLE
   IPOPT_LAPACK_FUNC(sgetrf, SGETRF)(&M, &N, a, &LDA, ipiv, &INFO);
#else
   IPOPT_LAPACK_FUNC(dgetrf, DGETRF)(&M, &N, a, &LDA, ipiv, &INFO);
#endif

   info = INFO;
#else

   std::string msg =
      "Ipopt has been compiled without LAPACK routine DGETRF, but options are chosen that require this dependency.  Abort.";
   THROW_EXCEPTION(LAPACK_NOT_INCLUDED, msg);
#endif

}

void IpLapackGetrs(
   Index         ndim,
   Index         nrhs,
   const Number* a,
   Index         lda,
   Index*        ipiv,
   Number*       b,
   Index         ldb
)
{
#ifdef IPOPT_HAS_LAPACK
   ipfint N = ndim, NRHS = nrhs, LDA = lda, LDB = ldb, INFO;
   char trans = 'N';

#ifdef IPOPT_SINGLE
   IPOPT_LAPACK_FUNC(sgetrs, SGETRS)(&trans, &N, &NRHS, a, &LDA, ipiv, b, &LDB,
                                    &INFO, 1);
#else
   IPOPT_LAPACK_FUNC(dgetrs, DGETRS)(&trans, &N, &NRHS, a, &LDA, ipiv, b, &LDB,
                                    &INFO, 1);
#endif

   DBG_ASSERT(INFO == 0);
#else

   std::string msg =
      "Ipopt has been compiled without LAPACK routine DGETRS, but options are chosen that require this dependency.  Abort.";
   THROW_EXCEPTION(LAPACK_NOT_INCLUDED, msg);
#endif

}

void IpLapackPpsv(
   Index         ndim,
   Index         nrhs,
   const Number* a,
   Number*       b,
   Index         ldb,
   Index&        info
)
{
#ifdef IPOPT_HAS_LAPACK
   ipfint N = ndim, NRHS = nrhs, LDB = ldb, INFO;
   char uplo = 'U';

#ifdef IPOPT_SINGLE
   IPOPT_LAPACK_FUNC(sppsv, SPPSV)(&uplo, &N, &NRHS, a, b, &LDB, &INFO);
#else
   IPOPT_LAPACK_FUNC(dppsv, DPPSV)(&uplo, &N, &NRHS, a, b, &LDB, &INFO);
#endif

   info = INFO;
#else

   std::string msg =
      "Ipopt has been compiled without LAPACK routine DPPSV, but options are chosen that require this dependency.  Abort.";
   THROW_EXCEPTION(LAPACK_NOT_INCLUDED, msg);
#endif
}

} // namespace Ipopt
#endif  /* ifdef IPOPT_HAS_LAPACK */
