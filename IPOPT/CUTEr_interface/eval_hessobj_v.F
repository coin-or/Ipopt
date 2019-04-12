C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
      subroutine EVAL_HESSOBJ_V(TASK, N, X, M, VIN, VOUT, DAT, IDAT)
C
C    $Id: eval_hessobj_v.F 531 2004-03-11 01:31:07Z andreasw $
C
C     TASK = 0: reevaluate Hessian (because X or LAM changed)
C            1: do not need to reevaluate - but X and LAM are still
C               set to the correct values
      implicit none
      integer TASK, N, M
      double precision  VIN(N), X(N), VOUT(N)
      double precision DAT(*)
      integer IDAT(*)

#include "CUTE.INC"
      double precision lam(CUTE_MMAX)

      call DCOPY(CUTE_M, 0d0, 0, lam, 1)

      call EVAL_HESSLAG_V(0, N, X, M, lam, VIN, VOUT, DAT, IDAT)

      return
      end
