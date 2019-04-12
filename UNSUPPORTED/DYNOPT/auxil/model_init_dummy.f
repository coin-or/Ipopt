C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine MODEL_INIT(NZ, NY, NU, NP, T, Z, DZ, Y, U, P, F,
     1     libname)
      implicit none
      double precision T
      integer NZ, NY, NU, NP
      double precision Z(NZ), DZ(NZ), Y(NY), U(NU), P(NP), F(NZ+NY)
      character*(*) libname
C
C     This is a dummy routine for get_start.f in case ADOL-C is not used
C     and no tape has to be written
C
      return
      end
