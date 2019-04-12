C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine MODEL( NZ, NY, NU, NP, T, Z, ZDOT, Y, U, P, F )
C
C     Model of Roscoe's Tank
C     Andreas Waechter  05-25-00
C                       08-25-01  Fortran model version
C                       08-28-01  ADOL-C version
C
      implicit none

      integer NZ, NY, NU, NP
      double precision T        ! Time
      double precision Z(NZ)    ! Differential variables
      double precision ZDOT(NZ) ! Derivatives of differential variables
      double precision Y(NY)    ! Algebraic variables
      double precision U(NU)    ! Control Variables
      double precision P(NP)    ! Parameters (might be subject to disturbances)
      double precision F(NZ+NY) ! RHS od DAE system
C
C     names for variables
C
      double precision h, hdot, hbar, pres, p_atm, F0, F1, rho, g
      double precision Area, CV

      integer l_z, l_y, l_u, l_p, l_ae, l_de

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
      l_z  = l_z + 1
      h    = Z   (l_z)
      hdot = ZDOT(l_z)
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C                            Algebraic Variables
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      l_y  = l_y + 1
      F1   = Y(l_y)
C
      l_y  = l_y + 1
      pres = Y(l_y) * 0.1d6
C
      l_y  = l_y + 1
      hbar = Y(l_y)
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C                             Control Variables
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      l_u  = l_u + 1
      CV   = U(l_u)
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
      F(l_de)  = (F0-F1)/Area - hdot
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C                             Algebraic Equations
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      l_ae     = l_ae + 1
      F(l_ae)  = F1**2 - (CV**2)*(pres-p_atm)
C
      l_ae     = l_ae + 1
      F(l_ae)  = (pres - p_atm - rho*g*h) / 0.1d6
C
      l_ae     = l_ae + 1
      F(l_ae)  = h - hbar
C
      return
      end
