C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine MODEL_HV( NRHS, NZ, NY, NU, NP, T, Z, ZDOT, Y, U,
     1     P, LAM, RHS, RES, LDRS)
C
C     Model of Roscoe's Tank
C     Andreas Waechter  08-25-01  Fortran model version
C
      implicit none

      integer NRHS, NZ, NY, NU, NP
      double precision T        ! Time
      double precision Z(NZ)    ! Differential variables
      double precision ZDOT(NZ) ! Derivatives of differential variables
      double precision Y(NY)    ! Algebraic variables
      double precision U(NU)    ! Control Variables
      double precision P(NP)    ! Parameters
      double precision LAM(NZ+NY) ! weights for individual Hessians
      integer LDRS              ! leading dimension of RHS, RES
      double precision RHS(LDRS, NRHS) ! vectors to be multiplied
      double precision RES(LDRS, NRHS) ! result
C
C     names for variables
C
      double precision h, hdot, hbar, pres, p_atm, F0, F1, rho, g
      double precision Area, CV

      integer l_z, l_y, l_u, l_p, l_ae, l_de
      integer l_h, l_hdot, l_F1, l_pres, l_hbar, l_CV

      integer j

C =============================================================================
C
C                           Constants in model
C
C =============================================================================

      include 'TANK_CONSTS.INC'

C =============================================================================
C
C                      Initialize Product Result to Zero
C
C =============================================================================

      do j = 1, NRHS
         call DCOPY(2*NZ+NY+NU+NP, 0.d0, 0, RES(1,j), 1)
      enddo

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
      l_ z   = l_z + 1
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
C     Loop over all RHS
C

      do j = 1, NRHS
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
C         F(l_de)  = (F0-F1)/Area - hdot
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C                             Algebraic Equations
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
         l_ae     = l_ae + 1
C         F(l_ae)  = F1**2 - (CV**2)*(pres-p_atm)
         RES(l_F1  ,j) = RES(l_F1  ,j) +
     1        LAM(l_ae)*( 2.d0             )*RHS(l_F1  ,j)
         RES(l_CV  ,j) = RES(l_CV  ,j) +
     1        LAM(l_ae)*(-2.d0*(pres-p_atm))*RHS(l_CV  ,j)
         RES(l_CV  ,j) = RES(l_CV  ,j) +
     1        LAM(l_ae)*(-2.d0*CV          )*RHS(l_pres,j) * 0.1d6
         RES(l_pres,j) = RES(l_pres,j) +
     1        LAM(l_ae)*(-2.d0*CV          )*RHS(l_CV  ,j) * 0.1d6
C
         l_ae     = l_ae + 1
C         F(l_ae)  = pres - p_atm - rho*g*h
C
         l_ae     = l_ae + 1
C         F(l_ae)  = h - hbar

C
C     End Loop
C
      enddo
C
      return
      end
