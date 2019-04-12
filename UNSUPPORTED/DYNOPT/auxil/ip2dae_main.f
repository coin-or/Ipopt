C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      program IP2NLP_MAIN
C
C     Main program for Dynamic Optimization using IPOPT and elemental
C     decompositon
C
C     $Id: ip2dae_main.f 531 2004-03-11 01:31:07Z andreasw $
C
C     Author:    Andreas Waechter
C                Department of Chemical Engineering
C                Carnegie Mellon University
C                Pittsburgh, PA
C                USA
C
C     Date       10-11-00
C
      implicit none
C
      include 'DYNAUX.INC'
C
      integer IW(LIW_IPOPT)
      double precision RW(LRW_IPOPT)
D      integer PIW, PRW, MYALLOC

      integer NX, M, NLB, NUB
      integer ILB(NXMAX), IUB(NXMAX)
      double precision X(NXMAX), BNDS_L(NXMAX), BNDS_U(NXMAX)
      double precision V_L(NXMAX), V_U(NXMAX), LAM(NXMAX)
      double precision times, timef

      integer nargs
      double precision args(50)
      character*20 cargs(50)

      integer i, ierr, iter, iele, icol, lx, j
      logical ex
C
C     Determine size of problem (in COMMON block IP2NLP), data about
C     discretization, bounds, initial conditions, starting point, and if
C     necessary initialize model
C
      call GET_START(NX, X, NLB, ILB, BNDS_L, NUB, IUB, BNDS_U, ierr)
      if( ierr.ne.0 ) then
         write(*,*) 'GET_START returns ierr = ',ierr
         stop
      endif
C
C     compute number of equality constraints
C
      M = NZ + NELE*(NCOL*(NZ+NY) + NZ + NAC)
cC
cC     Add bounds to obtain SS in last element...
cC
c      inquire(file = 'SS_BOUNDS', exist = ex)
c      if( ex ) then
c         write(*,*) 'Adding bounds to obtain SS in last element...'
c         lx = NZ + (NELE-1)*(NCOL*(NZ+NY)+NZ)
c         do i = 1, NCOL
c            do j = 1, NZ-1
c               NLB         = NLB + 1
c               ILB(NLB)    = lx + j
c               BNDS_L(NLB) = -1d-2
c               NUB         = NUB + 1
c               IUB(NUB)    = lx + j
c               BNDS_U(NUB) =  1d-2
c            enddo
c            lx = lx + NY + NZ
c         enddo
c      endif
C
C     If the file 'X.start' exists, read initial point from there
C
      inquire(file='X.start', exist = ex)
      if( ex ) then
         open(10,file='X.start')
         write(*,*) 'Reading initial point from X.start'
         do i = 1, NX
            read(10,'(d23.15)') X(i)
         enddo
         close(10)
      endif
C
C     Set numerical parameters
C
      nargs = 1
      cargs(1) = 'ifull'
      args(1) = 0.d0
      call INITPARAMS(ierr, nargs, args, cargs)
C
C     Choose the first element as the one based on  which parition and
C     factorization are to be determined
C
      CRIT_ELE = 1
C
      write(*,*) 'N = ',NX,' M = ',M,' NLB = ',NLB,' NUB = ',NUB
C
C     call IPOPT to solve optimization problem
C
      call TIMER(times)
      call IPOPT(NX, X, M, NLB, ILB, BNDS_L, NUB, IUB, BNDS_U, V_L, V_U,
     2     LAM, LRW_IPOPT, RW, LIW_IPOPT, IW, iter, ierr)

C      PRW = MYALLOC(8*LRW_IPOPT)
C      PIW = MYALLOC(4*LIW_IPOPT)
C      call IPOPT(NX, X, M, NLB, ILB, BNDS_L, NUB, IUB, BNDS_U, V_L, V_U,
C     2     LAM, LRW_IPOPT, %VAL(PRW), LIW_IPOPT, %VAL(PIW), iter, ierr)

      call TIMER(timef)
      write(*,*) 'total CPU time in IPOPT = ',timef-times
      if( ierr.ne.0 ) then
         write(*,*) 'IPOPT returns ierr = ', ierr
c         stop
      endif
C
C     For now, just write solution into file
C
      lx = 0
      open(7,file='X.sol',status='unknown')
      do i = 1, NZ
         write(7,3000) X(lx+i),' Z',i,0
      enddo
      lx = lx + NZ
      do iele = 1, NELE
         do icol = 1, NCOL
            do i = 1, NZ
               write(7,3000) X(lx+i),'dZ',i,iele,icol
            enddo
            lx = lx + NZ
            do i = 1, NY
               write(7,3000) X(lx+i),' Y',i,iele,icol
            enddo
            lx = lx + NY
         enddo
         do i = 1, NZ
            write(7,3000) X(lx+i),' Z',i,iele
         enddo
         lx = lx + NZ
      enddo
      do iele = 1, NELE
         do icol = 1, NCOL
            do i = 1, NU
               write(7,3000) X(lx+i),' U',i,iele,icol
            enddo
            lx = lx + NU
         enddo
      enddo
      do i = 1, NP
         write(7,3000) X(lx+i),' P',i
      enddo

 3000 format(d23.15,'  ',a2,'(',i3,')  ',i3,'  ',i2)

c      do i = 1, NX
c         write(7,'(d23.15,i8)') X(i),i
c      enddo
      close(7)
C
C     That's it?
C
      end
