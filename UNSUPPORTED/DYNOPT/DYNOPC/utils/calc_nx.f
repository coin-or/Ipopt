C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: calc_nx.f 531 2004-03-11 01:31:07Z andreasw $
      integer function CALC_NX(NZ, NY, NU, NP, NELE, NCOL)
!DEC$ ATTRIBUTES DLLEXPORT :: CALC_NX
C
C     Simply function to calculate size of discretized X vector for given
C     dimensions
C
C     Authors:    Yidong Lang, Andreas Waechter   09-29-01
C
      implicit none
      integer NZ, NY, NU, NP, NELE, NCOL

      CALC_NX = NZ + NELE*(NCOL*(NZ+NY+NU) + NZ) + NP
      return
      end
