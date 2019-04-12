C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: const_prof.f 531 2004-03-11 01:31:07Z andreasw $
      subroutine CONST_PROF(Z0, DZ, ZB, Y0, YB, U0, UB, P0, PB, X,
     1     IERR)
C
C     Copy initial estimates from *.cmp file into X to obtain
C     constant initial profile.
C
C     Author:  Andreas Waechter     10-07-01
C
      implicit none

      include 'DYNAUX.INC'
      include 'DYNOPC.INC'

      double precision Z0(*), DZ(*), ZB(2,*)
      double precision Y0(*), YB(2,*)
      double precision U0(*), UB(2,*)
      double precision P0(*), PB(2,*)
      double precision X(*)
      integer IERR

      integer lx, i, j, CALC_NX

      IERR = 0

C     Check dimensions

      lx = CALC_NX(NZ,NY,NU,NP,NELE,NCOL)
      if( lx.gt.NXMAX ) then
         IERR = -1
         write(*,*) 'Error on const_prof: NXMAX too small.'
         return
      endif
C
C     Copy steady state into X to obtain starting point for IPOPT
C
      lx  = 0

      call DCOPY(NZ, Z0, 1, X(lx+1), 1)
      lx = lx + NZ

      do i = 1, NELE
         do j = 1, NCOL
            call DCOPY(NZ, DZ, 1, X(lx+1), 1) ! zdot
            lx = lx + NZ
            call DCOPY(NY, Y0, 1, X(lx+1), 1) ! y
            lx = lx + NY
         enddo
         call DCOPY(NZ, Z0, 1, X(lx+1), 1) ! z
         lx = lx + NZ
      enddo

      do i = 1, NELE
         do j = 1, NCOL
            call DCOPY(NU, U0, 1, X(lx+1), 1) ! u
            lx = lx + NU
         enddo
      enddo

      call DCOPY(NP, P0, 1, X(lx+1), 1) ! p

      return
      end
