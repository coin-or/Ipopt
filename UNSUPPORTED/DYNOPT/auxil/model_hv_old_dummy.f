C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine MODEL_HV_OLD(NZ, NY, NU, NP, T, Z,
     1     DZ, Y, U, P, LAM, VZ, VDZ, VY, VU, VP,
     2     HVZ, HVDZ, HVY, HVU, HVP)
      implicit none
      double precision T
      integer NZ, NY, NU, NP
      double precision Z(NZ), DZ(NZ), Y(NY), U(NU), P(NP), LAM(NZ+NY)
      double precision VZ(NZ), VDZ(NZ), VY(NY), VU(NU), VP(NP)
      double precision HVZ(NZ), HVDZ(NZ), HVY(NY), HVU(NU), HVP(NP)
C
C     This is a dummy routine for lagvec.f in case ADOL-C is not used
C
      return
      end
