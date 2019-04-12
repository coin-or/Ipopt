C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine GET_START(NX, X, NLB, ILB, BNDS_L, NUB, IUB, BNDS_U,
     1     IERR)
C
C     Read start information from files.
C
C     Version for Fortran and ADOL-C interface
C
C     $Id: get_start.f 531 2004-03-11 01:31:07Z andreasw $
C
C     Author:  Andreas Waechter
C              Department of Chemical Engineering
C              Carnegie Mellon University
C              Pittsburgh, PA
C              USA
C
C     Date:    10-11-00
C              09-05-01    added constant controls and parameters
C
      implicit none
C
      include 'DYNAUX.INC'
C
      integer NX
      double precision X(NXMAX)
      integer NLB
      integer ILB(NXMAX)
      double precision BNDS_L(NXMAX)
      integer NUB
      integer IUB(NXMAX)
      double precision BNDS_U(NXMAX)
      integer IERR
C
      double precision Z(NZMAX), Z_L(NZMAX), Z_U(NZMAX), DZ(NZMAX)
      double precision Y(NYMAX), Y_L(NYMAX), Y_U(NYMAX)
      double precision U(NUMAX), U_L(NUMAX), U_U(NUMAX)
      double precision P(NPMAX), P_L(NPMAX), P_U(NPMAX)
      double precision F(NZMAX+NYMAX)
      logical SS(NZMAX)

      integer i, j, k, lx, lp, lu, ndegu, NU_OPT, NP_OPT
      double precision FinalTime, delta
      logical ex
      character*80 line

      IERR = 0
C
C     Read number of elements and collocation points
C
      open(10, file='SIZE.DAT',status='old')
      read(10,8010) NELE
      if( NELE.gt.NELEMAX ) goto 9000
      read(10,8010) NCOL
      if( NCOL.gt.NCOLMAX ) goto 9000
      read(10,8010) ndegu
      close(10)
C
C     Read steady state
C
      open(10,file='SS.DAT',status='old')
      read(10,8010) NZ
      if( NZ.gt.NZMAX ) goto 9000
      read(10,8010) NY
      if( NY.gt.NYMAX ) goto 9000
      read(10,8010) NU
      if( NU.gt.NUMAX ) goto 9000
      read(10,8010) NP
      if( NP.gt.NPMAX ) goto 9000
      if( NZ.gt.0 ) read(10,8000)
      do i = 1, NZ
         read(10,8020) Z(i), Z_L(i), Z_U(i)
      enddo
      if( NY.gt.0 ) read(10,8000)
      do i = 1, NY
         read(10,8020) Y(i), Y_L(i), Y_U(i)
      enddo
      if( NU.gt.0 ) read(10,8000)
      do i = 1, NU
         read(10,8020) U(i), U_L(i), U_U(i)
      enddo
      if( NP.gt.0 ) read(10,8000)
      do i = 1, NP
         read(10,8020) P(i), P_L(i), P_U(i)
      enddo
C
C     Read indices and values of fixed U's
C
      do i = 1, NU
         IU_PROF(i) = 0
      enddo
      read(10,8000)
      read(10,8010) NU_PROF
      do i = 1, NU_PROF
         read(10,8010) j
         if( j.gt.NU ) goto 9000
         IU_PROF(j) = 1
         read (10,8020) U_PROF(j,1)
         call DCOPY(NELE, U_PROF(j,1), 0, U_PROF(j,2), 1)
      enddo
C
C     Read indices and values of fixed P's
C
      do i = 1, NP
         IP_FIX(i) = 0
      enddo
      read(10,8000)
      read(10,8010) NP_FIX
      do i = 1, NP_FIX
         read(10,8010) j
         if( j.gt.NP ) goto 9000
         IP_FIX(j) = 1
         read (10,8020) P_FIX(j)
      enddo

      close(10)

      call DCOPY(NZ, 0.d0, 0, DZ, 1)

      NU_OPT = NU - NU_PROF
      NP_OPT = NP - NP_FIX

      NX = NZ + NELE*(NCOL*(NZ+NY+NU_OPT) + NZ) + NP_OPT
      if( NX.gt.NXMAX ) goto 9000
C
C     Initialize additional constraints
C
      call ADDCON_INIT(NZ, NY, NU_OPT, NCOL, ndegu, NAC)
C
C     Read finite element bounds
C
      open(10,file='TI.DAT',status='old')
      read(10,8000) line
      if( line.eq.'all' ) then
         do i = 1, NELE+1
            read(10,8020) TI(i)
         enddo
         close(10)
      elseif( line.eq.'linear' ) then
         read(10,8020) TI(1)
         read(10,8020) TI(NELE+1)
         do i = 2, NELE
            TI(i) = TI(1) + dble(i-1)/dble(NELE)*(TI(NELE+1)-TI(1))
         enddo
      elseif( line.eq.'airsep' ) then
         read(10,8020) FinalTime
         write(*,*) 'FinalTime = ',FinalTime
         delta = 2*FinalTime/(NELE*(NELE+1))
         do i=0,NELE
            ti(i+1)=delta*(i*(i+1)/2)
         enddo
      else
         write(*,*) 'Wrong key word in TI.DAT'
         IERR = -8
         goto 9999
      endif
C
C     Read initial conditions
C
      open(10,file='ZINIT.DAT',status='old')
      do i = 1, NZ
         read(10,8020) ZINIT(i)
      enddo
      close(10)
C
C     Create tape for ADOL-C  (dummy call if ADOL-C is not used)
C
      call MODEL_INIT(NZ, NY, NU, NP, TI(1), Z, DZ, Y, U, P, F,
     1     './model')
C
C     Copy steady state into X to obtain starting point for IPOPT
C
      lx  = 0
      NLB = 0
      NUB = 0

      call DCOPY(NZ, Z, 1, X(lx+1), 1)
      lx = lx + NZ

      do i = 1, NELE
         do j = 1, NCOL
            call DCOPY(NZ, DZ, 1, X(lx+1), 1) ! zdot
            lx = lx + NZ
            call DCOPY(NY, Y, 1, X(lx+1), 1) ! y
            do k = 1, NY
               if( Y_L(k).gt.-VLARGE ) then
                  NLB = NLB + 1
                  ILB(NLB) = lx + k
                  BNDS_L(NLB) = Y_L(k)
               endif
               if( Y_U(k).lt. VLARGE ) then
                  NUB = NUB + 1
                  IUB(NUB) = lx + k
                  BNDS_U(NUB) = Y_U(k)
               endif
            enddo
            lx = lx + NY
         enddo
         call DCOPY(NZ, Z, 1, X(lx+1), 1) ! z
         do k = 1, NZ
            if( Z_L(k).gt.-VLARGE ) then
               NLB = NLB + 1
               ILB(NLB) = lx + k
               BNDS_L(NLB) = Z_L(k)
            endif
            if( Z_U(k).lt. VLARGE ) then
               NUB = NUB + 1
               IUB(NUB) = lx + k
               BNDS_U(NUB) = Z_U(k)
            endif
         enddo
         lx = lx + NZ
      enddo

      do i = 1, NELE
         do j = 1, NCOL
            lu = 0
            do k = 1, NU
               if( IU_PROF(k).eq.0 ) then
                  lu = lu + 1
                  X(lx+lu) = U(k)
                  if( U_L(k).gt.-VLARGE ) then
                     NLB = NLB + 1
                     ILB(NLB) = lx + lu
                     BNDS_L(NLB) = U_L(k)
                  endif
                  if( U_U(k).lt. VLARGE ) then
                     NUB = NUB + 1
                     IUB(NUB) = lx + lu
                     BNDS_U(NUB) = U_U(k)
                  endif
               endif
            enddo
            lx = lx + NU_OPT
         enddo
      enddo

      lp = 0
      do k = 1, NP
         if( IP_FIX(k).eq.0 ) then
            lp = lp + 1
            X(lx+lp) = P(k)
            if( P_L(k).gt.-VLARGE ) then
               NLB = NLB + 1
               ILB(NLB) = lx + lp
               BNDS_L(NLB) = P_L(k)
            endif
            if( P_U(k).lt. VLARGE ) then
               NUB = NUB + 1
               IUB(NUB) = lx + lp
               BNDS_U(NUB) = P_U(k)
            endif
         endif
      enddo
      lx = lx + NP_OPT
C
C     Add bounds to obtain SS in last element...
C
      inquire(file = 'SS_BOUNDS', exist = ex)
      if( ex ) then
         write(*,*) 'Adding bounds to obtain SS in last element...'
         open(10,file='SS_BOUNDS', status='old')
         do i = 1, NZ
            SS(i) = .false.
         enddo
 200     read(10,8010,end=210) i
            if( i.gt.NZ ) goto 9000
            SS(i) = .true.
         goto 200
 210     close(10)
         lx = NZ + (NELE-1)*(NCOL*(NZ+NY)+NZ)
         do i = 1, NCOL
            do j = 1, NZ
               if( SS(j) ) then
                  NLB         = NLB + 1
                  ILB(NLB)    = lx + j
CORIG                  BNDS_L(NLB) = -1d-4*dmax1(1d0,dabs(ZINIT(j)))
                  BNDS_L(NLB) = -1d-3*dmax1(1d0,dabs(ZINIT(j)))
                  NUB         = NUB + 1
                  IUB(NUB)    = lx + j
CORIG                  BNDS_U(NUB) =  1d-4*dmax1(1d0,dabs(ZINIT(j)))
                  BNDS_U(NUB) =  1d-3*dmax1(1d0,dabs(ZINIT(j)))
               endif
            enddo
            lx = lx + NY + NZ
         enddo
      endif
C
C     That should be it...
C
      goto 9999

 8000 format(a)
 8010 format(i10)
 8020 format(3d23.15)

 9000 continue
      write(*,*) 'some maximal dimension exceeded - abort.'
      IERR = -1
      goto 9999

 9999 continue
      return
      end
