// Copyright (C) 2005, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Andreas Waechter              IBM    2005-12-25

#include "IpoptConfig.h"
#include "IpLapack.hpp"
#include "IpTypes.h"

#ifdef IPOPT_SINGLE
#define IPOPT_LAPACK_FUNCP(name,NAME) IPOPT_LAPACK_FUNC(s ## name,S ## NAME)
#else
#define IPOPT_LAPACK_FUNCP(name,NAME) IPOPT_LAPACK_FUNC(d ## name,D ## NAME)
#endif

#ifdef IPOPT_HAS_LAPACK
// Prototypes for the LAPACK routines
extern "C"
{
   /** LAPACK Fortran subroutine XPOTRS. */
   void IPOPT_LAPACK_FUNCP(potrs, POTRS)(
      char*           uplo,
      ipindex*        n,
      ipindex*        nrhs,
      const ipnumber* A,
      ipindex*        ldA,
      ipnumber*       B,
      ipindex*        ldB,
      ipindex*        info,
      int             uplo_len
   );

   /** LAPACK Fortran subroutine XPOTRF. */
   void IPOPT_LAPACK_FUNCP(potrf, POTRF)(
      char*     uplo,
      ipindex*  n,
      ipnumber* A,
      ipindex*  ldA,
      ipindex*  info,
      int       uplo_len
   );

   /** LAPACK Fortran subroutine XSYEV */
   void IPOPT_LAPACK_FUNCP(syev, SYEV)(
      char*     jobz,
      char*     uplo,
      ipindex*  n,
      ipnumber* A,
      ipindex*  ldA,
      ipnumber* W,
      ipnumber* WORK,
      ipindex*  LWORK,
      ipindex*  info,
      int       jobz_len,
      int       uplo_len
   );

   /** LAPACK Fortran subroutine XGETRF. */
   void IPOPT_LAPACK_FUNCP(getrf, GETRF)(
      ipindex*      m,
      ipindex*      n,
      ipnumber*     A,
      ipindex*      ldA,
      ipindex*      IPIV,
      ipindex*      info
   );

   /** LAPACK Fortran subroutine XGETRS. */
   void IPOPT_LAPACK_FUNCP(getrs, GETRS)(
      char*           trans,
      ipindex*        n,
      ipindex*        nrhs,
      const ipnumber* A,
      ipindex*        ldA,
      ipindex*        IPIV,
      ipnumber*       B,
      ipindex*        ldB,
      ipindex*        info,
      int             trans_len
   );

   /** LAPACK Fortran subroutine XPPSV. */
   void IPOPT_LAPACK_FUNCP(ppsv, PPSV)(
      char*           uplo,
      ipindex*        n,
      ipindex*        nrhs,
      const ipnumber* A,
      ipnumber*       B,
      ipindex*        ldB,
      ipindex*        info
   );
}

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
   ipindex N = ndim, NRHS = nrhs, LDA = lda, LDB = ldb, INFO;
   char uplo = 'L';

   IPOPT_LAPACK_FUNCP(potrs, POTRS)(&uplo, &N, &NRHS, a, &LDA, b, &LDB, &INFO, 1);
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
   ipindex N = ndim, LDA = lda, INFO;

   char UPLO = 'L';

   IPOPT_LAPACK_FUNCP(potrf, POTRF)(&UPLO, &N, a, &LDA, &INFO, 1);

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
   ipindex N = ndim, LDA = lda, INFO;

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
   ipindex LWORK = -1;
   Number WORK_PROBE;
   IPOPT_LAPACK_FUNCP(syev, SYEV)(&JOBZ, &UPLO, &N, a, &LDA, w,
                                  &WORK_PROBE, &LWORK, &INFO, 1, 1);
   DBG_ASSERT(INFO == 0);

   LWORK = (ipindex) WORK_PROBE;
   DBG_ASSERT(LWORK > 0);

   Number* WORK = new Number[LWORK];
   for (Index i = 0; i < LWORK; i++)
   {
      WORK[i] = i;
   }
   IPOPT_LAPACK_FUNCP(syev, SYEV)(&JOBZ, &UPLO, &N, a, &LDA, w,
                                  WORK, &LWORK, &INFO, 1, 1);

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
   ipindex M = ndim, N = ndim, LDA = lda, INFO;

   IPOPT_LAPACK_FUNCP(getrf, GETRF)(&M, &N, a, &LDA, ipiv, &INFO);

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
   ipindex N = ndim, NRHS = nrhs, LDA = lda, LDB = ldb, INFO;
   char trans = 'N';

   IPOPT_LAPACK_FUNCP(getrs, GETRS)(&trans, &N, &NRHS, a, &LDA, ipiv, b, &LDB,
                                    &INFO, 1);

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
   ipindex N = ndim, NRHS = nrhs, LDB = ldb, INFO;
   char uplo = 'U';

   IPOPT_LAPACK_FUNCP(ppsv, PPSV)(&uplo, &N, &NRHS, a, b, &LDB, &INFO);

   info = INFO;
#else

   std::string msg =
      "Ipopt has been compiled without LAPACK routine DPPSV, but options are chosen that require this dependency.  Abort.";
   THROW_EXCEPTION(LAPACK_NOT_INCLUDED, msg);
#endif
}

} // namespace Ipopt
#endif  /* ifdef IPOPT_HAS_LAPACK */
