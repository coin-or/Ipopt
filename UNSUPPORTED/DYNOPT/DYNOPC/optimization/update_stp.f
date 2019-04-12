C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: update_stp.f 531 2004-03-11 01:31:07Z andreasw $
      subroutine UPDATE_STP(X, TI_STP, NELE_STP, NCOL_STP, IERR)
C
C     Check whether starting point from file fits the current discretization,
C     otherwise update the starting point my intrapolation.
C
C     Authors:  Yidong Lang, Andreas Waechter   10-02-01
C
      implicit none
C
      double precision X(*)
      double precision TI_STP(*)
      integer NELE_STP, NCOL_STP, IERR
C
      include 'DYNAUX.INC'
      include 'DYNOPC.INC'
      include 'DAE2NLP.INC'
!DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/, /DYNOPC/, /DAENLP/
C
      integer i, nx, nx_old, lx, lu, iele, icol
      double precision x_old(NXMAX), dt, t
      double precision zdummy(NZMAX), ydummy(NYMAX), udummy(NUMAX)

      integer CALC_NX

      IERR = 0
C     Check, if intrapolation has to be done

      if( NELE_STP.ne.NELE ) goto 100
      if( NCOL_STP.ne.NCOL ) goto 100
      do i = 1, NELE+1
         if( TI_STP(i).ne.TI(i) ) goto 100
      enddo

C     Data is Ok, return
      return

C     Perform intrapolation
 100  continue

      write(*,105)
 105  format(/,6x,'Discretization of time horizon has changed',/,
     1     /,6x,'---> obtain new starting point by intrapolation.'/)

      nx = CALC_NX(NZ, NY, NU, NP, NELE, NCOL)
      if( nx.gt.NXMAX) then
         write(*,110) nx,NXMAX
 110     format(/,6x,'Error in UPDATE_STP: NX = ',i6,
     1        ', but NXMAX = ',i6)
         IERR = -1
         return
      endif

      nx_old = CALC_NX(NZ, NY, NU, NP, NELE_STP, NCOL_STP)

C     Copy old starting point
      call DCOPY(nx_old, X, 1, x_old, 1)

C     initialize data in COMMON block /DAENLP/ containing collocation points
C     'RHO(i)'

      call INITD2N(NCOL)

C     First NZ entries are Ok (initial conditions)
      lx = NZ
      lu = NZ + NELE*(NCOL*(NZ+NY)+NZ)

C     Loop over all finite elements and collocation points to construct X
      do iele = 1, NELE
         dt = TI(iele+1) - TI(iele)
         do icol = 1, NCOL
            t = TI(iele) + RHO(icol)*dt
            call APPSLN(t, NZ, NY, NU, NELE_STP, NCOL_STP, TI_STP,
     1           x_old, zdummy, X(lx+1), X(lx+NZ+1), X(lu+1))
            lx = lx + NZ + NY
            lu = lu + NU
         enddo
         call APPSLN(TI(iele+1), NZ, NY, NU, NELE_STP, NCOL_STP, TI_STP,
     1           x_old, X(lx+1), zdummy, ydummy, udummy)
         lx = lx + NZ
      enddo

C     Finally the unchanged parameters
      call DCOPY(NP, x_old(nx_old-NP+1), 1, X(lu+1), 1)

C     That's it!

      return
      end
