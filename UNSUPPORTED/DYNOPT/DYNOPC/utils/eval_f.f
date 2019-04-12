C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: eval_f.f 531 2004-03-11 01:31:07Z andreasw $
      subroutine EVAL_F(N, X, F)
C
C     Objective function for IPOPT
C     Linear terms that we defined in DYNOPC.INC
C
C     Authors:  Yidong Lang, Andreas Waechter     09-29-01
C
      implicit none
      integer N
      double precision F, X(N)

      include 'DYNAUX.INC'
      include 'DYNOPC.INC'
!DEC$ ATTRIBUTES DLLEXPORT :: /DYNAUX/, /DYNOPC/

      integer i, iex, j

      F = 0.d0
C
C     Linear terms (see DYNOPC.INC)
C
      do i = 1, NE_IN_OBJ
         iex = (E_INDEX(i)-1)*(NCOL*(NZ+NY) + NZ)
         do j = 1, NZ_IN_OBJ
            F = F + AIJZ(j,i)*X(iex + Z_INDEX(j))
         end do
CTODO: NEED TO WORK ON THIS what for E_INDEX = 1????
         iex = iex - NY
         do j = 1, NY_IN_OBJ
            F = F + BIJY(j,i)*X(iex+Y_INDEX(j))
         end do
      end do
      iex = NZ + NELE*(NCOL*(NZ+NY)+NZ)+NELE*NCOL*NU 
      do i = 1, NP_IN_OBJ
         F = F + CIP(i)*X(iex+P_INDEX(i))
      end do

      return
      end
