C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
      subroutine EVAL_HESSCON_V(TASK, N, X, M, LAM, VIN, VOUT,
     1     DAT, IDAT)
C
C    $Id: eval_hesscon_v.F 531 2004-03-11 01:31:07Z andreasw $
C
C     TASK = 0: reevaluate Hessian (because X or LAM changed)
C            1: do not need to reevaluate - but X and LAM are still
C               set to the correct values
      implicit none
      integer TASK, N, M
      double precision LAM(M), VIN(N), X(N), VOUT(N)
      double precision DAT(*)
      integer IDAT(*)

#include "CUTE.INC"
      double precision lam2(CUTE_MMAX)
      double precision tmp(CUTE_NMAX)

      if( M.gt.CUTE_MMAX .or. N.gt.CUTE_NMAX) then
         write(*,*) 'N or M too large in eval_hesscon_v.'
         stop
      endif

      call EVAL_HESSLAG_V(0, N, X, M, LAM, VIN, VOUT, DAT, IDAT)

      call DCOPY(CUTE_M, 0d0, 0, lam2, 1)
      call EVAL_HESSLAG_V(0, N, X, M, lam2, VIN, tmp, DAT, IDAT)
      call DAXPY(CUTE_N, -1d0, tmp, 1, VOUT, 1)

      return
      end
