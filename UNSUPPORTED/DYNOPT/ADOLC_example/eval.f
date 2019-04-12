C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine EVAL_F(N, X, F)

      implicit none
      integer n
      double precision f, x(N)

      include 'DYNAUX.INC'

      F = X(NZ + NELE*(NCOL*(NZ+NY)+NZ))

      return
      end

      subroutine EVAL_G(N, X, G)
C
C    $Id: eval.f 531 2004-03-11 01:31:07Z andreasw $

      implicit none
      integer n
      double precision g(n), x(N)
      include 'DYNAUX.INC'

      call DCOPY(N, 0.d0, 0, G, 1)
      G(NZ + NELE*(NCOL*(NZ+NY)+NZ)) = 1d0

      return
      end

      subroutine EVAL_HESSOBJ_V(TASK, N, X, M, VIN, VOUT)

      implicit none
      integer TASK, N, M
      double precision  VIN(N), X(N), VOUT(N)

      call DCOPY(N, 0d0, 0, VOUT, 1)

      return
      end
