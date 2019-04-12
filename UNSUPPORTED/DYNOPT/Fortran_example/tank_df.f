C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine MODEL_DF( NZ, NY, NU, NP, T, Z, ZDOT, Y, U, P,
     1     DF, LDDF )
C
C     Model of Roscoe's Tank - Jacobian
C     Andreas Waechter  08-25-01  Fortran model version
C
      implicit none

      integer NZ, NY, NU, NP
      double precision T        ! Time
      double precision Z(NZ)    ! Differential variables
      double precision ZDOT(NZ) ! Derivatives of differential variables
      double precision Y(NY)    ! Algebraic variables
      double precision U(NU)    ! Control Variables
      double precision P(NP)    ! Parameters
      integer          LDDF     ! Leading dimension of DF
      double precision DF(LDDF,2*NZ+NY+NU+NP) 
                                ! Jacobian of DAE system
C
C     names for variables
C
      double precision h, hdot, hbar, pres, p_atm, F0, F1, rho, g
      double precision Area, CV

      integer l_z, l_y, l_u, l_p, l_ae, l_de
      integer l_h, l_hdot, l_F1, l_pres, l_hbar, l_CV

C =============================================================================
C
C                           Constants in model
C
C =============================================================================

      include 'TANK_CONSTS.INC'

C =============================================================================
C
C                  Copy Z, Y, U, P into variables with names
C
C =============================================================================

C
C     Initialize counters for variables
C
      l_z  = 0
      l_y  = 0
      l_p  = 0
      l_u  = 0
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C                           Differential Variables
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      l_z    = l_z + 1
      h      = Z   (l_z)
      l_h    = l_z
      hdot   = ZDOT(l_z)
      l_hdot = NZ+l_z
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C                            Algebraic Variables
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      l_y    = l_y + 1
      F1     = Y(l_y)
      l_F1   = 2*NZ+l_y
C
      l_y    = l_y + 1
      pres   = Y(l_y) * 0.1d6
      l_pres = 2*NZ+l_y
C
      l_y    = l_y + 1
      hbar   = Y(l_y)
      l_hbar = 2*NZ+l_y
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C                             Control Variables
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      l_u    = l_u + 1
      CV     = U(l_u)
      l_CV   = 2*NZ+NY+l_u
C
C =============================================================================
C
C                                  M O D E L
C
C =============================================================================
C
C
C     Initialize counters for equations
C
      l_de = 0
      l_ae = NZ
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C                            Differential Equations
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      l_de     = l_de + 1
C      F(l_de)  = (F0-F1)/Area - hdot
      DF(l_de,l_F1  ) = -1.d0/Area
      DF(l_de,l_hdot) = -1.d0

C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C                             Algebraic Equations
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      l_ae     = l_ae + 1
C      F(l_ae)  = F1**2 - (CV**2)*(pres-p_atm)
      DF(l_ae,l_F1  ) = 2.d0*F1
      DF(l_ae,l_CV  ) = -2.d0*CV*(pres-p_atm)
      DF(l_ae,l_pres) = -(CV**2) * 0.1d6
C
      l_ae     = l_ae + 1
C      F(l_ae)  = (pres - p_atm - rho*g*h) / 0.1d6
      DF(l_ae,l_pres) = 1.d0 ! / 0.1d6
      DF(l_ae,l_h   ) = -rho*g / 0.1d6
C
      l_ae     = l_ae + 1
C      F(l_ae)  = h - hbar
      DF(l_ae,l_h   ) = 1.d0
      DF(l_ae,l_hbar) = -1.d0
C
      return
      end
