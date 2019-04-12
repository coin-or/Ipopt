C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: eval_g.f 531 2004-03-11 01:31:07Z andreasw $
      subroutine EVAL_G(N, X, G)
C
C     Gradient of objective function for IPOPT
C     Linear terms that we defined in DYNOPC.INC
C
C     Authors:  Yidong Lang, Andreas Waechter     09-29-01
C
      implicit none
      integer N
      double precision G(N), X(N)

      include 'DYNAUX.INC'
      include 'DYNOPC.INC'
!DEC$ ATTRIBUTES DLLEXPORT :: /DYNAUX/, /DYNOPC/

      integer i, iex, j

      call DCOPY(N, 0.d0, 0, G, 1)

      do i = 1, NE_IN_OBJ
         iex = (E_INDEX(i)-1)*(NCOL*(NZ+NY) + NZ)
         do j = 1, NZ_IN_OBJ
            G(iex+Z_INDEX(j)) = AIJZ(j,i)
         end do
         iex = iex - NY
         do j = 1, NY_IN_OBJ
            G(iex+Y_INDEX(j)) = BIJY(j,i)
         end do
      end do

      iex = NZ + NELE*(NCOL*(NZ+NY)+NZ)+NELE*NCOL*NU
      do i = 1, NP_IN_OBJ
          G(iex+P_INDEX(i)) = CIP(i)
      end do

      return
      end
